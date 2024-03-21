# Optimization_NelderMead
A Matlab project for optimization of Rosenbrock function (3 variables) using Nelder Mead method with visualization

## Equation:
f(x_1,x_2,x_3 )= (1 - x_1 )^2  + (2 -x_2 )^2+ 65*((x_2-(x_1)^2))^2+ 65*((x_3-(x_2)^2))^2

### Start point: [1 -2 5]

#### Convergence:
〖(∑_(i=1)^(n+1)▒(f_i-f_c )^2/(n+1))〗^(1/2)≤ε 

stop algorithm when ε=1e-6

### result:

Number of Iterations
reflection    152
contraction   29
shrinkage	    0
expansion	    35
Total         216

                    Optimum
          X	                               F
     [1.365, 1.866, 3.484]               0.1516
                 	



