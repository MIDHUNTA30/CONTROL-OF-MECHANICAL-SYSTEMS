% MATLAB program for Sliding Mode Control of PMSM 
% main function
function main
clear all;
close all;
% System parameters
Rs=3;Ld=.01;Lq=.01;
J=.001;B=.001;P=4;f=.5;Kt=.1;

xk=[0,0];
zk=[0,0];
uk=0; sk=0;
T=0.001;
for k=1:1:20000
t(k)=k*T;
x1r(k)=5*sin(.1*k*T)+25;
time=[0 T];
[tt,xx]=ode45(@(t,x) pmsm(t,x,uk),time,xk);
xk=xx(length(xx),:); 
x1(k)=xk(1);
x2(k)=xk(2);
s(k)=x1r(k)-x1(k);
[tt,zz]=ode45(@(t,z) smo(t,z,sk),time,zk);  % Estimating s dot
zk=zz(length(zz),:); 
z1(k)=zk(1);
z2(k)=zk(2);

w(k)=10000*sign(z1(k))+1800*sign(z2(k));
v(k)=(J/P*Kt)*((B/J)*((B/J)*x1(k)+(P*Kt/J)*x2(k))+w(k));
u(k)=Lq*((Rs/Lq)*x2(k)+(f/Lq)*x1(k)+v(k));

uk=u(k);
sk=s(k);
end

figure(1);
subplot(211);
plot(t,x1,'r',t,x1r,'k','linewidth',2);
xlabel('Time(sec)', 'Interpreter','latex');ylabel('$\omega$', 'Interpreter','latex');
legend('$\omega$','$\omega_{r}$', 'Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(212);
plot(t,s,'k',t,z1,'r','linewidth',2);
xlabel('Time(sec)', 'Interpreter','latex');ylabel('$s$', 'Interpreter','latex');
legend('$s$','Estimated $s$', 'Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg f1
figure(2);
subplot(211);
plot(t,z2,'r','linewidth',2)
xlabel('Time(sec)',  'Interpreter','latex');ylabel('Estimated $\dot{s}$','Interpreter','latex');
%legend('Estimated  $\dot{s}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(212);
plot(t,u,'r','linewidth',2)
xlabel('Time(sec)',  'Interpreter','latex');ylabel('$u$', 'Interpreter','latex');
%legend('u');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg f2

% state equation: PMSM
function dx=pmsm(t,x,u)
dx=zeros(2,1);
Rs=3;Ld=.01;Lq=.01;
J=.001;B=.001*(1+.001*t);P=4;f=.5;Kt=.1;Tl=0*sin(t);

dx(1)=-(B/J)*x(1)+(P*Kt/J)*x(2)-(P/J)*Tl;
dx(2)=-(Rs/Lq)*x(2)-(f/Lq)*x(1)+(1/Lq)*u;

% state equation: sliding mode observer
function dz=smo(t,z,s)
dz=zeros(2,1);
dz(1)=-10*((z(1)-s)^.5)*sign(z(1)-s)+z(2);
dz(2)=-50*sign(z(1)-s);
