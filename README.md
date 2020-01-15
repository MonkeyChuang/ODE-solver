# ODE Solver

## Introduction
``ode_solver.m`` provides several methods to solve system of nonlinear differential equations of the form dy/dt=f(t,y).

## Arguments
Its input and output formats is 'similar' (but not the same) to that of `ode45` in MATLAB.
### INPUT
1. `method`:
    - Format: string
     - Description: Specify the method used to solve the problem. Currently contains 
        - "forwardEuler" -- forward Euler method
        - "backwardEuler" -- backward Euler method
        - "AdamBashforth2" -- Adam-Bashforth 2-step method
        - "AdamMoulton2"-- Adam-Moulton 2-step method
        - "RungeKutta4" -- Fourth Runge-Kutta method
        - "RungeKutta45" -- Runge-Kutta-Fehlberg method (adaptive step size)
1. `f`: 
    - Format: function handle
    - Description: function to solve. **This function MUST return a COLUMN vector**
1. `tspan`:
    - Format: vector
    - Description: mesh points (or partition) of interval of itegration. It should be given in the form [t0,t1,...,tf]. **!!This argument is different to that of ode45!!**
1. `y_init`:
    - Format: vector
    - Description: initial conditions.

### OUTPUT
1. `t`:
   - Format: column vector
   - Description: the same as input 'tspan'
1. `y`:
   - Format: array
   - Description: Solution. Each row y(i,:) represents the value of solution at time t(i).  