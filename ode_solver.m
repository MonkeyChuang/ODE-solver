% ==========================Documentation==========================
% 'ode_solver' provides several methods to solve system of nonlinear 
% differential equations of the form dy/dt=f(t,y).
% 
% Its input and output formats is 'similar' (but not the same) 
% to that of 'ode45':
%   *Input arguments
%       'method':
%           Format: string
%           Description: Specify the method used to solve the problem.
%                       Currently contains 
%                       "forwardEuler" -- forward Euler method
%                       "backwardEuler" -- backward Euler method
%                       "AdamBashforth2" -- Adam-Bashforth 2-step method
%                       "AdamMoulton2"-- Adam-Moulton 2-step method
%                       "RungeKutta4" -- Fourth Runge-Kutta method
%                       "RungeKutta45" -- Runge-Kutta-Fehlberg method (adaptive step size)
%       'f': 
%           Format: function handle
%           Description: function to solve. 
%                       This function MUST return a COLUMN vector 
%       'tspan':
%           Format: vector
%           Description: mesh points (or partition) of interval of itegration.
%                       It should be given in the form [t0,t1,...,tf].
%                       !!This argument is different to that of ode45!!
%       'y_init':
%           Format: vector
%           Description: initial conditions.
%   *Output arguments
%       't':
%           Format: column vector
%           Description: the same as input 'tspan'
%       'y':
%           Format: array
%           Description: Solution. 
%                       Each row y(i,:) represents the value of solution at time t(i).  
% ===================================================================
function [t,y] = ode_solver(method,f,tspan,y_init)
    if method=="forwardEuler"
        [t,y]=forwardEuler(f,tspan,y_init);
    elseif method=="backwardEuler"
        [t,y]=backwardEuler(f,tspan,y_init);
    elseif method=="AdamBashforth2"
        [t,y]=AdamBashforth2(f,tspan,y_init);
    elseif method=="AdamMoulton2"
        [t,y]=AdamMoulton2(f,tspan,y_init);
    elseif method=="RungeKutta4"
        [t,y] = RK4(f,tspan,y_init);
    elseif method=="RungeKutta45"
        [t,y] = RK45(f,tspan,y_init);
    end
end
function [t,y] = forwardEuler(f,tspan,y_init)
    t = tspan;
    n = length(y_init);                                 % n: dimension of domain & range of f.
    N = length(t);                                      % N: mesh size(including t0)
    
    y = zeros(N,n);                                     % y: array that stores solution at each t
    y(1,:) = y_init';       
    for k=2:N
        dt = t(k)-t(k-1);                               % We do not assume that mesh points are equaldistant
        y(k,:) = y(k-1,:)+f(t(k-1),y(k-1,:))'*dt;
    end
end
function [t,y] = backwardEuler(f,tspan,y_init)
    t = tspan;
    n = length(y_init);                                 % n: dimension of domain & range of f.
    N = length(t);                                      % N: mesh size(including t0)
    
    y = zeros(N,n);                                     % y: array that stores solution at each t
    y(1,:) = y_init';
    
    for k=2:N
        dt = t(k)-t(k-1);                               % We do not assume that mesh points are equaldistant
        G = @(u) u-y(k-1,:)'-f(t(k),u)*dt;                 % input/output of G: column vector
        %DG = @(u) eye(n)-Df(t(k),u)*dt;
        y_predict = y(k-1,:)' + f(t(k-1),y(k-1,:))*dt; % use 1 step forward Euler to predict y(k,:)
        
        %[y(k,:),~,iterations] = Newton(G,DG,y_predict,1000,1e-08);
        [y(k,:),~,iterations] = Broyden(G,y_predict,1000,1e-08);
    end
end
function [t,y] = AdamBashforth2(f,tspan,y_init)
    t = tspan;
    N = length(t);                                      % N: mesh size(including t0)
    n = length(y_init);                                 % n: dimension of domain & range of f.
    
    y = zeros(N,n);                                     % y: array that stores solution at each t
    dt = t(2)-t(1);
    y(1,:) = y_init';
    y(2,:) = y(1,:)+f(t(1),y(1,:))'*dt;                 % initial guess of y2 using forward Euler
    for k=3:N
        dt = t(k)-t(k-1);                               % We do not assume that mesh points are equaldistant
        A = 3*f(t(k-1),y(k-1,:))/2;
        B = -f(t(k-2),y(k-2,:))/2;
        y(k,:) = y(k-1,:)+(A+B)'*dt;
    end    
end
function [t,y] = AdamMoulton2(f,tspan,y_init)
    t = tspan;
    N = length(t);                                      % N: mesh size(including t0)
    n = length(y_init);                                 % n: dimension of domain & range of f.
    
    y = zeros(N,n);                                     % y: array that stores solution at each t
    dt = t(2)-t(1);
    y(1,:) = y_init';
    y(2,:) = y(1,:)+f(t(1),y(1,:))'*dt;                 % initial guess of y2 using forward Euler
    for k=3:N
        dt = t(k)-t(k-1);                               % We do not assume that mesh points are equaldistant
        A = 2*f(t(k-1),y(k-1,:))/3;
        B = -f(t(k-2),y(k-2,:))/12;
        G = @(u) u-y(k-1,:)'-( 5*f(t(k),u)/12+A+B )*dt;
        y_predict = y(k-1,:)' + f(t(k-1),y(k-1,:))*dt;
        
        [y(k,:),~,iterations] = Broyden(G,y_predict,1000,1e-08);
    end    

end
function [x,converged,iterations] = Broyden(f,x0,maxiter,tol)
    % =====
    % f and x0 must be column vector
    % output x is also a column vector
    % =====
    iter = 0;
    n = length(x0);      
    x_new = x0;
    
    % Approximate Df(x0) using forward difference with step size h.
    h=1e-4;
    Df= zeros(n,n);
    e = zeros(n,1);
    e(1) = 1;
    for i=1:n
        Df(:,i) = (f(x0+h*e)-f(x0))/h;
        e = circshift(e,1);
    end
    % Start iterative step
    while norm(f(x_new))>=tol && iter<maxiter
        x_old = x_new;
        % When n(dimension of Df) is large, one can use Sherman-Morrison
        % formula to compute Df^(-1)f(x_old) instead of Df\f(x_old).
        x_new = x_old-Df\f(x_old);
        Df = Df+ (f(x_new)-f(x_old)-Df*(x_new-x_old))*(x_new-x_old)'/norm(x_new-x_old)^2;
        iter = iter+1;
    end
    if norm(f(x_new))<tol
        converged = true;
    else
        converged = false;
    end
    iterations = iter;
    x = x_new;
end

function [t,y] = RK4(f,tspan,y_init)
    t = tspan;
    n = length(y_init);                                 % n: dimension of domain & range of f.
    N = length(t);                                      % N: mesh size(including t0)
    
    y = zeros(N,n);                                     % y: array that stores solution at each t
    y(1,:) = y_init';
    for k=2:N
        dt = t(k)-t(k-1);
        k1 = dt*f(t(k-1),y(k-1,:))';            % row vector
        k2 = dt*f(t(k-1)+dt/2,y(k-1,:)+k1/2)';  % row
        k3 = dt*f(t(k-1)+dt/2,y(k-1,:)+k2/2)';  % row
        k4 = dt*f(t(k),y(k-1,:)+k3)';           % row
        
        y(k,:) = y(k-1,:)+(k1 + 2*k2 + 2*k3 + k4)/6;
    end
end
% RK45 �u�ݪ��D�_�l�ɶ� a �P�פ�ɶ� b 
function [t,y] = RK45(f,tspan,y_init)
    a = tspan(1);                                       % a: initial time
    b = tspan(end);                                     % b: end time
    
    y(1,:) = y_init';                                   % y: array that stores solution at each t
    t(1) = a;
    
    tol = 1e0;
    dt = 3.9;
    
    max_iter = 5000;
    k = 1;
    while t(k)<b && k<max_iter
        dt = min(dt,b-t(k));
        
        w = y(k,:);                                 % For the sake of clarity, we change the notation
        T = t(k);                                   % Same as above.
        k1 = dt*f(T,w)';                                                        %row vector
        k2 = dt*f(T+dt/4 ,w+k1/4)';                                             %row
        k3 = dt*f(T+3*dt/8 ,w+3*k1/32+9*k2/32)';                                %row
        k4 = dt*f(T+12*dt/13,w+1932*k1/2197-7200*k2/2197+7296*k3/2197)';        %row
        k5 = dt*f(T+dt,w+439*k1/216-8*k2+3680*k3/513-845*k4/4104)';             %row
        k6 = dt*f(T+dt/2,w-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40)';   %row
        
        %w1 = w+25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        %w2 = w+16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        %R = abs(w1-w2)/dt;
        R = norm(k1/360-128*k3/4275-2197*k4/75240+k5/50+2*k6/55)/dt;
        delta = 0.84*(tol/R)^(1/4);
        
        if R<=tol
            %fprintf("Step %d (T=%f) -- Accept! R=%f, delta=%f, dt=%e\n ",k,T,R,delta,dt)
            fprintf("Step %d (T=%f) dt=%e\n ",k,T,dt)
            t(k+1) = t(k)+dt;
            y(k+1,:) = w+25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
            k = k+1;
        else
            %fprintf("Step %d (T=%f)-- Reject! R=%f, delta=%f, dt=%e\n ",k,T,R,delta,dt)
            fprintf("Step %d (T=%f) dt=%e\n ",k,T,dt)
        end
        % The if-else block prevent dt from being expanding by large delta
        
        if delta<=0.1
            dt = 0.1*dt;
        elseif delta>=4
            dt = 4*dt;
        else
            dt = delta*dt;
        end
        

    end
end

