clear 
close all

c = 9.8;
y0 = [pi/16;0];
dt = 0.125;
tspan = 0:dt:10;

f = @(t,y)[y(2) ; -c*sin(y(1))];

figure
[tspan,y] = ode_solver("RungeKutta4",f,tspan,y0);
%[~,y_exact] = ode45(f,tspan,y0);
%plot(tspan,y(:,1),tspan,y_exact(:,1),'r')
plot(tspan,y(:,1))



