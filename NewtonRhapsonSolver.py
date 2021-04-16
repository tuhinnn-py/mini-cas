# TODO: Use the Newton-Rhapson Method for finding out the real roots of a function f.
# Note: Newton-Rhapson Method has a second order of convergence (2)

# Method : b = a - f(a) / f'(a), where a is the initial approximation to the root
# Method : Replace a with b for the next iteration
# Note: Newton-Rhapson Method uses one initial approximation to the root

# Newton-Rhapson Method can be proved in one of the following two ways:
#   When the curve meets the X-axis, the curve can be approximated with a straight line (a1.x + a0) in an infinitesimally small neighborhood
#   f(x) = a1.x + a0 ==> f'(x) = a1 ==> a0 = f(x) - f'(x).x ==> x = -a0 /a1 ==> x_new = x - f(x) / f'(x)

#   If x* is the exact real root of a function f, an approximation x can be written as x = x* - h, without loss of generality
#   Using Taylor's expansion Theorem of f at x_0: f(x) = f(x_0) + (x - x_0).f'(x_0) + O(x), where O(x) stands for higher order functions which can be ignored for the purpose of the proof
#   f(x*) = f(x_0) + (x* - x_0)f'(x_0) = 0, since f(x*) = 0 ==> x* - x_0 = -f(x_0) / f'(x_0) ==> x* = x_0 - f(x_0) / f'(x_0)

# Note: In an iterative numerical method, we find the next approximation using a function g(x), where you could think of g as a generator function which generates the next(hopefully, more accurate) approximation of a root to a function f
# Note: In Newton-Rhapson Method, g(x) = x - f(x) / f'(x)
# Note: The two very important properties of the generator function are namely:
#           g(x*) = x*
#           g'(x) < tan(pi/4), In this case it's 0

# Note: Stop iteration in one of two situations:
#           |f(x)| < epsilon
#           |x(k + 1) - x(k)| < epsilon

# Note: Set manual approximations while calculating roots of f using init_approx
# Note: Instead of using manual approximations, set approximations of f using Bisection Method through self.set_approximations(*args, **kwargs)

from typing import Callable, Tuple, List, Sequence, Dict
from NumericalMethodSolver import NumericalMethodSolver

class NewtonRhapsonSolver(NumericalMethodSolver):
    def __init__(self, func : Callable, der_func : Callable, intervals : Sequence = None, epsilon : float = 1e-5) -> None:
        super(NewtonRhapsonSolver, self).__init__(func, intervals, epsilon)
        self.der_func = der_func

    @staticmethod
    def calculate_next_value(func : Callable, der_func : Callable, init_x : float, cache_value : bool = False) -> float:
        f_init_x = func(init_x)
        der_f_init_x = der_func(init_x)

        _ = f_init_x / der_f_init_x
        x_new = init_x - _
        
        if cache_value:
            # Provision for Chebysev's Method
            return x_new, f_init_x, der_f_init_x, _
        else:
            return x_new
    
    def set_roots(self, init_approx : float = None) -> None:
        if init_approx == None:
            assert not self.approximations == None, "Approximations need to be set before calculating roots by calling self.set_approximations(*args, **kwargs)"
            self.roots = []
            for approximation in self.approximations.keys():
                x_0 = approximation
                f_x_k = self.func(x_0)
                while(abs(f_x_k) >= self.epsilon):
                    x_new = self.calculate_next_value(self.func, self.der_func, x_0)
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
                    x_new = self.calculate_next_value(self.func, self.der_func, x_0)
                    x_0 = x_new
                    f_x_k = self.func(x_new)

                self.roots.append(x_0)
        
