rng(0);

%% Simulación
ti = 0; %tiempo inicial 
tf = 360; %tiempo final

v0 = [1;2;0.5]; %condiciones inciales tanques
p0 = [1/3;1/3;1/3]; %condiciones iniciales RD
x0 = [v0;p0];

tspan = [ti,tf];

beta = 0.01;
[t,x]=ode45(@(t,x)ODEx(t,x,beta),tspan,x0);

beta = 0.1;
[t2,x2]=ode45(@(t,x)ODEx(t,x,beta),tspan,x0);


%% Resultados
fig = figure;
fig.Position(3:4) = [290 300];
hold on
plot(t/60,x(:,1:3),'LineWidth',2)
plot(t2/60,x2(:,1),'LineStyle',':','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(t2/60,x2(:,2),'LineStyle',':','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(t2/60,x2(:,3),'LineStyle',':','Color',[0.9290 0.6940 0.1250],'LineWidth',2)
xlabel('Tiempo (minutos)','Interpreter','latex','FontSize',12)
ylabel('Volumen de agua ($m^3$)','Interpreter','latex','FontSize',12)
legend('tanque fuente 1','tanque fuente 2','tanque receptor','Interpreter','latex','FontSize',12)
axis([0,6,0,15])
grid on

fig = figure;
fig.Position(3:4) = [290 300];
hold on
plot(t/60,x(:,4:6)*100,'LineWidth',2)
plot(t2/60,x2(:,4)*100,'LineStyle',':','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(t2/60,x2(:,5)*100,'LineStyle',':','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(t2/60,x2(:,6)*100,'LineStyle',':','Color',[0.9290 0.6940 0.1250],'LineWidth',2)
axis([0,6,-50,90])
xlabel('Tiempo (minutos)','Interpreter','latex','FontSize',12)
ylabel('Proporci\''on de jugadores (\%)','Interpreter','latex','FontSize',12)
legend('estrategia 1 (\% apertura v\''alvula 1)','estrategia 2 (\% apertura v\''alvula 2)','estrategia 3','Interpreter','latex','FontSize',12)
yticks([-20 -10 0 10 20 30 40 50 60 70 80 90])
yticklabels({'','','0','10','20','30','40','50','60','70','80','90'})
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = ODEx(t,x,beta)

% preliminares
v1 = x(1);
v2 = x(2);
v3 = x(3);
p1 = x(4);
p2 = x(5);
p3 = x(6);

% parámetros

r = [0.5 0.25]; %lluvia
a1 = 0.1; %parámetro
a2 = 0.1;
a3 = 0.05;

% RD
P = 1; %masa poblacional
f1 = v1; %función de pago 1
f2 = v2; %función de pago 2
f3 = v3; %función de pago 3
fb = (p1*f1+p2*f2+p3*f3)/P; %pago promedio

% dinámicas
dv1 = r(1)-a1*p1*v1;
dv2 = r(2)-a2*p2*v2;
dv3 = a1*p1*v1+a2*p2*v2-a3*v3;
dp1 = beta*p1*(f1-fb);
dp2 = beta*p2*(f2-fb);
dp3 = beta*p3*(f3-fb);

dx(1,1) = dv1;
dx(2,1) = dv2;
dx(3,1) = dv3;
dx(4,1) = dp1;
dx(5,1) = dp2;
dx(6,1) = dp3;

end
