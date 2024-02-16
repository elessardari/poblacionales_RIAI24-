import numpy as np
import matplotlib.pyplot as plt
import math

np.random.seed(777)

N = 20               # Number of players (populations)
N = np.maximum(N, 2) # At least two players are considered
T = 40               # Total number of time slots

# Set competing time slots guaranteeing that nk >= 2 for all k
cond = False
while not cond:
    competing_time_slots = np.random.choice([1, 0], (N, T))
    competing_time_slots[-1, :] = np.ones(T) - competing_time_slots[0,:]
    nks = competing_time_slots.sum(1)
    if(np.min(nks) >= 2):
        if(np.min(np.sum(competing_time_slots, 0))>0):
            cond = True

competing_time_slots_aux = np.zeros_like(competing_time_slots)
for k in range(N):
    counter = 1
    for i in range(T):
        if(competing_time_slots[k, i] == 1):
            competing_time_slots_aux[k, i] = counter
            counter = counter + 1

# Define total number of variables
n = np.sum(nks)
print('Dimensional space (n):', n)

# Create Ck matrices
Cks = []
for k in range(N):
    Ck = np.zeros((T, nks[k]))
    temp = competing_time_slots[k, :]
    i = 0
    for j in range(T):
        if(temp[j] == 1):
            Ck[j, i] = 1
            i = i + 1
    Cks.append(Ck)
C = np.hstack(Cks)

# Define pseudo-gradient elements
Jbar = np.random.uniform(2, 4, (T, 1))
D = np.diag((np.random.uniform(0, 1, T)))
sqrt_D = np.sqrt(D)
alpha = np.random.uniform(1, 10, (n, 1))
beta = np.random.uniform(0, 1, (n, 1))

Rks = []
for k in range(N):
    Rks.append(np.dot(sqrt_D, Cks[k]))
R = np.hstack(Rks)

S = np.zeros((n, n))
pos = 0
for k in range(N):
    Ck = Cks[k]
    Sk = np.dot(Ck.T, np.dot(D, Ck))
    nk = nks[k]
    S[pos : pos+nk, pos : pos+nk] = Sk
    pos = pos + nk
S =  S + np.dot(R.T, R)

# Define pseudo-gradient mapping
def f(x):
    return -np.dot(S, x) - np.dot(C.T, Jbar) - (alpha * x) - beta 

# Define demand inequality constraints (Ax - b >= 0)
A = np.zeros((T, n))
for k in range(N):
    for i in range(T):
        if(competing_time_slots[k, i] == 1):
            pos = np.sum(nks[:k])
            A[i, pos + competing_time_slots_aux[k, i] - 1] = 1

b = np.random.uniform(2, 2.5, (T, 1))   # Upper bounds on energy demands

# Define simplex
mks = np.random.uniform(3, 4, (N, 1))   # Players' energy demands

# Define payoff vector
def p(x, mu):
    return f(x) - np.dot(A.T, mu)

# Define auxiliary matrix for Smith dynamics
G = np.zeros((n, n))
pos = 0
for k in range(N):
    block = np.ones((nks[k], nks[k])) - np.eye(nks[k])
    G[pos : pos + nks[k], pos : pos + nks[k]] = block
    pos = pos + nks[k]

# SIMULATION
Tsim = 10000
dt = 0.02

# Initialization
x0 = np.zeros((n, 1))
pos = 0
for k in range(N):
    x0k = np.ones((nks[k])) * (mks[k, 0]/nks[k])
    x0[pos : pos + nks[k], 0] = x0k
    pos = pos + nks[k]

mu0 = np.zeros((T, 1))

x = np.copy(x0)
mu = np.copy(mu0)

# Logs
xs = np.zeros((n, Tsim))
mus = np.zeros((T, Tsim))
cs = np.zeros((T, Tsim))
mass_CV = np.zeros((N, Tsim))
xs_u = np.zeros((n, Tsim))

# Main loop
for t in range(Tsim):
    xs[:, t] = x.reshape(-1)
    mus[:, t] = mu.reshape(-1)
    cs[:, t] = (np.dot(A, x) - b).reshape(-1)
    xs_u[:, t] = (x.reshape(-1)) / np.sqrt(np.dot(x.T, x))

    # PDM
    g = np.dot(A, x) - b
    dmu = np.maximum(g, 0) - mu*np.maximum(-g, 0)
    mu = mu + dt*dmu

    # EDM
    payoffs = p(x, mu)
    delta_payoffs = payoffs - payoffs.reshape(-1)
    x_pos_terms = np.maximum(delta_payoffs, 0) * (x.reshape(-1))
    x_neg_terms = np.maximum(-delta_payoffs, 0) * x
    dx = np.sum((x_pos_terms - x_neg_terms)*G, 1).reshape(n, 1)
    x = x + dt*dx

# KPIs
GNE_norm = np.sqrt(np.dot((x - x0).T, x - x0))
error = xs - x
KPI = np.sqrt((error * error).sum(0))/GNE_norm[0,0]

x_u = x / np.sqrt(np.dot(x.T, x))

thetas = np.arccos(np.clip(np.dot(x_u.T, xs_u), -1, 1))

x1s = (KPI * np.cos(thetas)).reshape(-1)
x2s = (KPI * np.sin(thetas)).reshape(-1)


print('MAX dmu:', np.max(np.abs(dmu)))
print('MAX dx:', np.max(np.abs(dx)))

# Figures
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 11})

t = dt*(np.arange(Tsim) + (1/dt))

plt.figure(1, figsize=(5, 4))
plt.plot(t, KPI, linewidth=2)
plt.grid()
plt.xscale('log')
plt.xlabel('Tiempo ' + r'$(t)$')
plt.ylabel(r'$\parallel\mathbf{x}(t) - \mathbf{x^*}\parallel_2\slash\parallel\mathbf{x}(0) - \mathbf{x^*}\parallel_2$')
plt.tight_layout(pad=0.2)
plt.savefig('resultados_JP1.pdf', bbox_inches='tight', pad_inches=0.2)

plt.figure(2, figsize=(5, 4))
for q in range(T):
    plt.plot(t, cs[q, :], linewidth=2)
plt.grid()
plt.xscale('log')
plt.ylabel('Restricciones de energia')
plt.xlabel('Tiempo ' + r'$(t)$')
plt.tight_layout(pad=0.2)
plt.savefig('resultados_JP2.pdf', bbox_inches='tight', pad_inches=0.2)

'''
plt.figure(4, figsize=(5, 4))
for q in range(T):
    plt.plot(t, mus[q, :], linewidth=2)
plt.grid()
plt.xscale('log')
plt.ylabel('Dual variables ' + r'$\mathbf{\mu}(t)$')
plt.xlabel('Time ' + r'$(t)$')
plt.tight_layout(pad=0.2)

plt.figure(5, figsize=(5, 4))
plt.plot(x1s, x2s)
plt.xlim((-1, 1))
plt.ylim((-1, 1))
'''
plt.show()
