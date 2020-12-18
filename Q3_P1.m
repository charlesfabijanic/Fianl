clc
clearvars
close all

Nx=10;
dx=15/Nx;


a=1;
dt=0.01;
T=zeros(1,Nx);

S=@(x)(-(x^2-4*x+2)*exp(-x));
T(end)=S(15);
time_itr=10000;
T_old=inf*ones(length(T));
k=0;
while norm(T_old-T)>10^-6
    x=0;
    k=k+1;
    for i=2:Nx-1
        x=x+dx;
        T_new(i)=T(i)+2*dt*((a/dx)*(T(i-1)-2*T(i)+T(i+1))+S(x));
    end
    T_old=T;
    T(1:Nx-1)=T_new;
   
end

plot(linspace(0,15,length(T)),T)
title('Nx=10')