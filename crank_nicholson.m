% Crank-Nicolson for the Heat Equation

%% Setup
close all;
N = 100;
L = 15;
alpha = 1;
dx = L/(N-1);
X = linspace(0,L,N)';

% initial condition
T = 0*X;
T(end)=0;
%% Spatial derivative operator
A = gallery('tridiag',N-2,1,-2,1);

%% inhomogeneous term
f = @(t) (-(X(2:N-1).^2-4*X(2:N-1)+2).*exp(-X(2:N-1)));

%% Time advancement
% Must solve the system 
%
% $$\left(I-\frac{\alpha dt}{2dx^2}A\right) T^{n+1} = \left(I+\frac{\alpha dt}{2dx^2}A\right) T^{n} + \frac{dt}{2}\left(f(t^{n+1})+f(t^n)\right)$$
% 

g = @(t,dt,T) (speye(N-2)-alpha*dt/2/(dx^2)*A)\(T(2:N-1)+ ...
              alpha*dt/2/(dx^2)*A*T(2:N-1) + dt*(f(t)+f(t+dt))/2);

%% Stable run
dt = 0.01;           % time step
t_final = 1000.0;        % final time
time = 0:dt:t_final;  % time array
pt = [0.0 0.5 1.0 1.5 t_final];   % desired plot times
pn = length(pt);          % number of desired plots
pc = 1;               % plot counter
rt = zeros(1,pn);
S = zeros(N,pn);      % solution storage
T_old=inf*ones(length(T));
k=0;
for t = time
    % plot storage
    k=k+1;
    if ( t >= pt(pc) )
        S(:,pc) = T;
        rt(pc) = t;
        pc = pc + 1;
        if (norm(T-T_old)<10^-3)
            break
        end

    end
    T_old=T;
    % time advancement
    T(2:N-1) = g(t,dt,T);
    if (norm(T-T_old)<10^-6)
        break
    end
end
norm(T-T_old)
%% Plot stable run
figure(1)
plot(X,T,'LineWidth',1)
xlabel('x','FontSize',14)
ylabel('T(x)','FontSize',14)
title('Nx=20')

