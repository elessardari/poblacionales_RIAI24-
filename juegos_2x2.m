rng(0);
%%
% Caso Dilema del Prisionero
% T>R>P>S

T = 0;
R = -1;
P = -5;
S = -20;

eq = (P-S)/((P-S)+(R-T));


J = [[R S];[T P]];
Tspan = 0:0.1:10;
p0 = rand(1,20);

for i=1:length(p0)
    [t,p] = ode23s(@(t,p) RD_ODE(t,p,J), Tspan, p0(i));
    
pDP(:,i) = p;
tDP(:,i) = t;
end

%%
% Caso Anti-coordinación
% T>R, S>P

% Este sería el caso de un HD game
% con v=2 y c=5
T = 0;
R = (2-5)/2;
P = 1;
S = 2;

eq = (P-S)/((P-S)+(R-T));

J = [[R S];[T P]];
Tspan = 0:0.1:10;
p0 = rand(1,20);

for i=1:length(p0)
    [t,p] = ode23s(@(t,p) RD_ODE(t,p,J), Tspan, p0(i));
    pAC(:,i) = p;
    tAC(:,i) = t;
end


%%
% Caso Coordinación
% R>T, P>S

% Este caso es el del 
% stag hunt porque además se
% tiene que T>P
T = 4;
R = 5;
P = 3;
S = 2;

eq = (P-S)/((P-S)+(R-T));

J = [[R S];[T P]];
Tspan = 0:0.1:10;
p0 = rand(1,20);

for i=1:length(p0)
    [t,p] = ode23s(@(t,p) RD_ODE(t,p,J), Tspan, p0(i));
    pC(:,i) = p;
tC(:,i) = t;
end


%%
% Caso Harmonía
% R>T, S>P

% En este caso, se tiene 
% R>T>S>P

T = 4;
R = 5;
P = 1;
S = 2;

eq = (P-S)/((P-S)+(R-T));


J = [[R S];[T P]];
Tspan = 0:0.1:10;
p0 = rand(1,20);

for i=1:length(p0)
    [t,p] = ode23s(@(t,p) RD_ODE(t,p,J), Tspan, p0(i));
    pH(:,i)  = p;
tH(:,i) = t;
    
end


fig = figure;
fig.Position(3:4) = [360 320];
plot(tDP,pDP, 'LineWidth',2)
%title('Replicator Dynamics para el Dilema del Prisionero', 'Interpreter', 'Latex', 'FontSize', 14)
xlabel('Tiempo (s)', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('$x_1(t)$', 'Interpreter', 'Latex', 'FontSize', 14)
grid on


fig = figure;
fig.Position(3:4) = [360 320];
plot(tAC,pAC, 'LineWidth',2)
%title("Replicator Dynamics para Anti-coordinaci\'on", 'Interpreter', 'Latex', 'FontSize', 14)
xlabel('Tiempo (s)', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('$x_1(t)$', 'Interpreter', 'Latex', 'FontSize', 14)
grid on


fig = figure;
fig.Position(3:4) = [360 320];
plot(tC,pC, 'LineWidth',2)
%title("Replicator Dynamics para Coordinaci\'on", 'Interpreter', 'Latex', 'FontSize', 14)
xlabel('Tiempo (s)', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('$x_1(t)$', 'Interpreter', 'Latex', 'FontSize', 14)
grid on


fig = figure;
fig.Position(3:4) = [360 320];
plot(tH,pH, 'LineWidth',2)
%title("Replicator Dynamics para Harmon\'ia", 'Interpreter', 'Latex', 'FontSize', 14)
xlabel('Tiempo (s)', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('$x_1(t)$', 'Interpreter', 'Latex', 'FontSize', 14)
grid on


 
function dx = RD_ODE(t,p,J)
    R = J(1,1);
    S = J(1,2);
    T = J(2,1);
    P = J(2,2);
    dx = p*(1-p)*(p*((P-S)+(R-T))-(P-S));
end