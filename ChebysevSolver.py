# TODO: Use the Chebysev's Method for finding out the real roots of a function f.
# Note: Chebysev's Method has a third order of convergence (3)

# Method : b = a - f(a) / f'(a) - (1/2).f(a)^2.f''(a)/f'(a)^3, where a is the initial approximation to the root
# Method : Replace a with b for the next iteration
# Note: Chebysev Method uses one initial approximation to the root

# Chebysev's Method can be proved in one of the following two ways:
#   When the function meets the X-axis, the function can be approximated with a quadratic curve(parabola) in an infinitesimally small neighborhood (a2.x^2 + a1.x + a0)
#   f(x_0) = a2.x_0^2 + a1.x_0 + a0                     -(0)
#   f'(x_0) = 2.a2.x_0 + a1 ==> f''(x_0) = 2.a2
#   a2 = f''(x_0) / 2                                   -(1)
#   a1 = f'(x_0) - f''(x_0).x_0                         -(2)
#   a0 = f(x_0) + (f''(x_0).x_0^2)/2 - f'(x_0).x_0      -(3)
#   Substituting the values of a2, a1, a0 from (1), (2) and (3) in (0) and simplifying it,
#   f(x_0) + (x - x_0)f'(x_0) + ((x - x_0)^2.f''(x_0))/2 = 0
#   Since, x = x* should satisfy this equation,
#   f(x_0) + (x* - x_0)f'(x_0) + ((x* - x_0)^2.f''(x_0))/2 = 0
#   x* - x_0 = -f(x_0)/f'(x_0) - (1/2).(x* - x_0)^2.f''(x_0)/f'(x_0)
#   Using the Newton-Rhapson Method for the latter (x* - x_0), we replace (x* - x_0) with -f(x_0)/f'(x_0)
#   x* - x_0 = -f(x_0)/f'(x_0) - (1/2).(f(x_0)/f'(x_0))^2.f''(x_0)/f'(x_0)
#   x* = x_0 - f(x_0)/f'(x_0) - (1/2).f(x_0)^2.f''(x_0)/f'(x_0)^3


#   Alternatively,
#   If x* is the exact real root of a function f, an approximation x_0 can be written as x_0 = x* - h, without loss of generality
#   Using Taylor's expansion Theorem of f at x_0: f(x) = f(x_0) + (x - x_0).f'(x_0) + ((x - x_0)^2.f''(x_0))/2 + O(x), where O(x) stands for higher order functions which can be ignored for the purpose of the proof
#   f(x*) = f(x_0) + (x* - x_0)f'(x_0) + ((x* - x_0)^2.f''(x_0))/2 = 0, since, f(x*) = 0
#   x* - x_0 = -f(x_0)/f'(x_0) - (1/2).(x* - x_0)^2.f''(x_0)/f'(x_0)
#   Using the Newton-Rhapson Method for the latter (x* - x_0), we replace (x* - x_0) with -f(x_0)/f'(x_0)
#   x* - x_0 = -f(x_0)/f'(x_0) - (1/2).(f(x_0)/f'(x_0))^2.f''(x_0)/f'(x_0)
#   x* = x_0 - f(x_0)/f'(x_0) - (1/2).f(x_0)^2.f''(x_0)/f'(x_0)^3

# Note: Stop iteration in one of two situations:
#           |f(x)| < epsilon
#           |x(k + 1) - x(k)| < epsilon

# Note: Set manual approximations while calculating roots of f using init_approx
# Note: Instead of using manual approximations, set approximations of f using Bisection Method through self.set_approximations(*args, **kwargs)

from typing import Callable, Tuple, List, Sequence, Dict
from NumericalMethodSolver import NumericalMethodSolver
from NewtonRhapsonSolver import NewtonRhapsonSolver

class ChebysevSolver(NumericalMethodSolver):
    def __init__(self, func : Callable, der_func : Callable, der_der_func : Callable, intervals : Sequence = None, epsilon : float = 1e-5) -> None:
        super(ChebysevSolver, self).__init__(func, intervals, epsilon)
        self.der_func = der_func
        self.der_der_func = der_der_func
    
    def set_roots(self, init_approx : float = None) -> None:
        if init_approx == None:
            assert not self.approximations == None, "Approximations need to be set before calculating roots by calling self.set_approximations(*args, **kwargs)"
            self.roots = []
            for approximation in self.approximations.keys():
                x_0 = approximation
                f_x_k = self.func(x_0)
                while(abs(f_x_k) >= self.epsilon):
                    x_new, f_init_x, der_f_init_x, _ = NewtonRhapsonSolver.calculate_next_value(self.func, self.der_func, x_0, True)
                    x_new -= 0.5 * f_init_x * f_init_x * self.der_der_func(x_0) / (der_f_init_x * der_f_init_x * der_f_init_x)
                    if (x_new - x_0) < self.epsilon:
                        x_0 = x_new
                        break
                    x_0 = x_new
                    f_x_k = self.func(x_new)

                self.roots.append(x_0)
        else:
            self.roots = []
            if isinstance(init_approx, float) or isinstance(init_approx, int):
                init_approx = [init_approx]
            for i in range(len(init_approx)):
                x_0 = init_approx[i]
                f_x_k = self.func(x_0)
                while(abs(f_x_k) >= self.epsilon):
                    x_new, f_init_x, der_f_init_x, _ = NewtonRhapsonSolver.calculate_next_value(self.func, self.der_func, x_0, True)
                    x_new -= 0.5 * f_init_x * f_init_x * self.der_der_func(x_0) / (der_f_init_x * der_f_init_x * der_f_init_x)
                    if (x_new - x_0) < self.epsilon:
                        x_0 = x_new
                        break
                    x_0 = x_new
                    f_x_k = self.func(x_new)

                self.roots.append(x_0)
        
