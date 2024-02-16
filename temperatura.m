rng(0);

%parámetros edificio
n = 4; %número de salones
c = 7.476*10^(4); %capacitancia térmica de cada salón
a = 123.6; %conductancia térmica de las paredes
P = 10000; %potencia de calefacción
Ts = [19; 20; 21; 22]; %referencias

%parámetros simulación
tf = 720; %tiempo final en minutos
tspan = [0,tf*60];
Ta = 1; %temperatura ambiente
T0 = Ta*ones(n,1); %temperaturas iniciales de los salones
x0 = P/(n+1)*ones(n+1,1);%condiciones iniciales RD (incluye zona ficticia)

z0 = [T0;x0];

%solución
[t,z] = ode45(@(t,z)TEMP(t,z,c,a,Ta,P,Ts,n),tspan,z0);

%resultados
T = z(:,1:n);
x = z(:,n+1:2*n+1);
xr = x(:,1:n);
xtot = sum(xr,2);

fig = figure;
fig.Position(3:4) = [290 260];
plot(t/60,T,'LineWidth',2)
legend('sal\''on 1: referencia 19$^{\circ}$C','sal\''on 2: referencia 20$^{\circ}$C','sal\''on 3: referencia 21$^{\circ}$C','sal\''on 4: referencia 22$^{\circ}$C','Interpreter','latex','FontSize',11)
xlabel('Tiempo (minutos)','Interpreter','latex','FontSize',13)
ylabel('Temperatura ($^{\circ}$C)','Interpreter','latex','FontSize',13)
yticks([5 10 15 19 20 21 22])
axis([0,tf,5,23])
grid

fig = figure;
fig.Position(3:4) = [290 260];
plot(t/60,xtot/1000,'LineWidth',2)
xlabel('Tiempo (minutos)','Interpreter','latex','FontSize',13)
ylabel('Potencia (kW)','Interpreter','latex','FontSize',13)
axis([0,tf,8,10])
grid

%dinámicas
function dz = TEMP(t,z,c,a,Ta,P,Ts,n)
    dst = 0;
    if t>200*60
        dst = -200;
    end

    Tamb = Ta;
    B = 20;
    beta = 0.0001;

    T1 = z(1);
    T2 = z(2);
    T3 = z(3);
    T4 = z(4);

    x1 = z(5);
    x2 = z(6);
    x3 = z(7);
    x4 = z(8);
    x5 = z(9);

    dT1 = a*((T2-T1) + (Ta-T1))/c + x1/c;
    dT2 = a*((T1-T2) + (T3-T2) + (Ta-T2))/c + x2/c;
    dT3 = a*((T2-T3) + (T4-T3) + (Ta-T3))/c + x3/c;
    dT4 = a*((T3-T4) + (Ta-T4))/c + x4/c + dst/c;

    f1 = Ts(1)-T1 + B;
    f2 = Ts(2)-T2 + B;
    f3 = Ts(3)-T3 + B;
    f4 = Ts(4)-T4 + B;
    f5 = B;

    fb = (x1*f1 + x2*f2 + x3*f3 + x4*f4 + x5*f5)/P;

    dx1 = beta*x1*(f1-fb);
    dx2 = beta*x2*(f2-fb);
    dx3 = beta*x3*(f3-fb);
    dx4 = beta*x4*(f4-fb);
    dx5 = beta*x5*(f5-fb);

    dz = [dT1;dT2;dT3;dT4;dx1;dx2;dx3;dx4;dx5];

end

