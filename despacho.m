rng(0);

%parámetros problema
n = 50; %número de generadores
G = 10000;
a = 0.01+0.004*rand(n,1);
b = 0.04*(1-2*rand(n,1));
c = 0.1+0.1*rand(n,1);


%parámetros simulación
tf = 30; %tiempo final en milisegundos
tspan = [0,tf/1000];
x0 = G/n*ones(n,1); %condiciones iniciales SD

%solución
[t,x] = ode45(@(t,x)Disp(t,x,n,a,b),tspan,x0);

%resultados
for i=1:n
    cost(:,i) = a(i)*x(:,i).^2 + b(i)*x(:,i) + c(i);
    mcost(:,i) = 2*a(i)*x(:,i) + b(i);
end
tcost = sum(cost,2);
xtot = sum(x,2);

fig = figure;
fig.Position(3:4) = [350 260];
plot(t*1000,tcost,'LineWidth',2)
xlabel('Tiempo (ms)','Interpreter','latex','FontSize',13)
ylabel('Costo total (\$)','Interpreter','latex','FontSize',13)
grid

fig = figure;
fig.Position(3:4) = [350 260];
plot(t*1000,xtot/1000,'LineWidth',2)
xlabel('Tiempo (ms)','Interpreter','latex','FontSize',13)
ylabel('Potencia generada (MW)','Interpreter','latex','FontSize',13)
axis([0,tf,0,20])
grid

fig = figure;
fig.Position(3:4) = [350 260];
plot(t*1000,mcost,'LineWidth',2)
xlabel('Tiempo (ms)','Interpreter','latex','FontSize',13)
ylabel('Costos marginales (\$/kW)','Interpreter','latex','FontSize',13)
grid

%dynamics
function dx = Disp(t,x,n,a,b)
    for i=1:n
        p(i) = -2*a(i)*x(i) - b(i);
    end

    for i=1:n
        dx(i,1) = 0;
        for j=1:n
        dx(i,1) = dx(i,1) + x(j)*max(p(i)-p(j),0) - x(i)*max(p(j)-p(i),0); 
        end
    end
end