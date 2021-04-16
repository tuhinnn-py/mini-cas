# TODO: Use the Multipoint Iteration Method for finding out the real roots of a function f.
# Note: The Multipoint Iteration Method has a third order of convergence (3)
# Note: |The Multipoint Iteration Method is used in cases where f'' does not exist but we would still like to have a third order of convergence

# The Multipoint Iteration Method can be proved in the following way:
# Using Taylor's expansion, f(x) = f(x_0) + (x - x_0).f'(x_0) + (x - x_0)^2.f''(x_0) + O(x), where O(x) stands for higher order functions which can be ignored for the purpose of this proof
# x - x_0 = -f(x_0) / (f'(x_0) + (1/2).(x - x_0).f''(x_0))                                    -(0)
# Expanding f'(x) about x_0 using Taylor's series about x_0,
# f'(x) = f'(x_0) + (x - x_0).f''(x_0) + Z(x), where Z(x) stands for higher order polynomials which can be ignored for the purpose of this proof
# f'(x_0 + (1/2).(x* - x_0)) = f'(x_0) + (1/2).(x* - x_0).f''(x_0) which is the denominator of (0)
# Therefore, x_new = x_0 - f(x_0)/f'(x_0 + (1/2).(x_new - x_0)), where x_new - x_0 can be approximated using the Newton-Rhapson Method

# Note: Stop iteration in one of two situations:
#           |f(x)| < epsilon
#           |x(k + 1) - x(k)| < epsilon

# Note: Set manual approximations while calculating roots of f using init_approx
# Note: Instead of using manual approximations, set approximations of f using Bisection Method through self.set_approximations(*args, **kwargs)

from typing import Callable, Tuple, List, Sequence, Dict
from NumericalMethodSolver import NumericalMethodSolver
import random

class MultipointSolver(NumericalMethodSolver):
    def __init__(self, func : Callable, der_func : Callable, intervals : Sequence = None, epsilon : float = 1e-5) -> None:
        super(MultipointSolver, self).__init__(func, intervals, epsilon)
        self.der_func = der_func
   
    def set_roots_1(self, init_approx : float = None) -> None:
        if init_approx == None:
            assert not self.approximations == None, "Approximations need to be set before calculating roots by calling self.set_approximations(*args, **kwargs)"
            self.roots = []
            for approximation in self.approximations.keys():
                x_0 = approximation
                f_x_k = self.func(x_0)
                while(abs(f_x_k) >= self.epsilon):
                    x_new = x_0 - f_x_k / self.der_func(x_0 - 0.5 * (f_x_k / self.der_func(x_0)))
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
                    x_new = x_0 - f_x_k / self.der_func(x_0 + 0.5 * (-f_x_k / self.der_func(x_0)))
                    x_0 = x_new
                    f_x_k = self.func(x_new)

                self.roots.append(x_0)
                
    def set_roots_2(self, init_approx : float = None) -> None:
        if init_approx == None:
            assert not self.approximations == None, "Approximations need to be set before calculating roots by calling self.set_approximations(*args, **kwargs)"
            self.roots = []
            for approximation in self.approximations.keys():
                x_0 = approximation
                f_x_k = self.func(x_0)
                while(abs(f_x_k) >= self.epsilon):
                    der_f_x_k = self.der_func(x_0)
                    x_inter = f_x_k / der_f_x_k
                    x_new = x_0 - x_inter - self.func(x_0 - x_inter) / der_f_x_k
                    
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
                    der_f_x_k = self.der_func(x_0)
                    x_inter = f_x_k / der_f_x_k
                    x_new = x_0 - x_inter - self.func(x_0 - x_inter) / der_f_x_k
                    
                    x_0 = x_new
                    f_x_k = self.func(x_new)

                self.roots.append(x_0)

    def set_roots(self, init_approx : float = None, use_method : int = 0) -> None:
        if use_method:
            self.set_roots_2(init_approx)
        else:
            self.set_roots_1(init_approx)
