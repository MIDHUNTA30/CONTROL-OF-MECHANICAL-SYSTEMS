% MATLAB program for Noninear MPC: Inverted pendulum system 
function main;
clear all;
close all
global N Np n m Ts x0
NT=500;N=10;n=4;m=1; Ts = 0.1;
xmin=-10;xmax=10;umin=-1;umax=1;
Q=eye(n); QN=5*Q; R=eye(m);
                       
x0=[0;0;1;0];
x=zeros(n,NT+1); x(:,1)=x0;
Xk=zeros(n*(N+1),1); Xk(1:n,1)=x0;
u=zeros(m,NT);
Uk=zeros(m*N,1);
zk=[Xk;Uk];

% constructing QX,RU,H
QX=Q;RU=R;
for i=1:N-1
  QX=blkdiag(QX,Q); RU=blkdiag(RU,R);
end
QX=blkdiag(QX,QN);
H=blkdiag(QX,RU);

% simulating system with MPC
for k=1:NT
   x0=x(:,k);  
   fun = @(z)z'*H*z;
   F=[];g=[];Feq=[];geq=[]; nonlcon=@nlcon;
   lb=[xmin*ones(1,(N+1)*n),umin*ones(1,N*m)];
   ub=[xmax*ones(1,(N+1)*n),umax*ones(1,N*m)];;
   z=fmincon(fun,zk,F,g,Feq,geq,lb,ub,nonlcon);
   u(:,k)=z((N+1)*n+1:(N+1)*n+m,1);
   x(:,k+1)=f(x(:,k),u(k));
   zk=z;
end    


% plotting response
figure(1)
time = (0:NT);
subplot(2,1,1)
plot(time,x(1,:),'r.-','LineWidth',.7)  
xlabel('$k$','Interpreter','latex');ylabel('${x}_{1}$','Interpreter','latex');
grid on
ax=gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,1,2)
plot(time,x(2,:),'k.-','LineWidth',.7)
xlabel('$k$','Interpreter','latex');ylabel('${x}_{2}$','Interpreter','latex');
grid on
ax=gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
figure(2)
time = (0:NT);
subplot(2,1,1)
plot(time,x(3,:),'r.-','LineWidth',.7)  
xlabel('$k$','Interpreter','latex');ylabel('${x}_{3}$','Interpreter','latex');
grid on
ax=gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,1,2)
plot(time,x(4,:),'k.-','LineWidth',.7)
xlabel('$k$','Interpreter','latex');ylabel('${x}_{4}$','Interpreter','latex');
grid on
ax=gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
figure(3)
stairs(time(1:end-1),u,'r.-','LineWidth',.7)
xlabel('$k$','Interpreter','latex');ylabel('${u}_{k}$','Interpreter','latex');
grid on
ax=gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
end

function xnext=f(x,u)
global NT N n m Ts x0     
g=9.8;M=2; mp=0.2; l=0.3;
xnext=x+Ts*[x(2);((M+mp)*g*sin(x(1)) + mp*l*cos(x(1))*sin(x(1))*(x(2)^2) + u*cos(x(1)))/(mp*l*(cos(x(1))^2)-(M+mp)*l);x(4); (mp*l*sin(x(1))*(x(2)^2)-mp*g*cos(x(1))*sin(x(1))+u)/(M+mp-mp*(cos(x(1))^2))];
end

function [c,ceq] =nlcon(z)
global NT N n m Ts x0 
c=[];
for i=1:N
    c1((i-1)*n+1:i*n,1)=z(i*n+1:(i+1)*n,1)-f(z((i-1)*n+1:i*n,1),z((N+1)*n+(i-1)*m+1:(N+1)*n+i*m));
end
ceq = [z(1:n)-x0;c1];
end


