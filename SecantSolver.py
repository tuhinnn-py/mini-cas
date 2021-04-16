# TODO: Use the Secant/Chord Method for finding out the real roots of a function f.
# Note: Secant Method has a superlinear order of convergence (1.618033), more commonly known as phi(golden ratio)

# Method : c = b - f(b) * (b - a) / (f(b) - f(a)), where a, b are intital approximations to the real root of the function
# Method : Replace a with b and b with c for the next iteration
# Note: Secant Method uses two initial approximations to the root

# Secant Method can be proved in the following way:
#   When the curve meets the X-axis, the curve can be approximated with a straight line (a1.x + a0) in an infinitesimally small neighborhood
#   Calculate the values of the coefficients(a0, a1) by solving the following system of equations : a1.a + a0 = f(a), a1.b + a0 = f(b): c = -a0 / a1

# Note: Stop iteration in one of two situations:
#           |f(x)| < epsilon
#           |x(k + 1) - x(k)| < epsilon

# Note: Set manual approximations while calculating roots of f using init_approx
# Note: Instead of using manual approximations, set approximations of f using Bisection Method through self.set_approximations(*args, **kwargs)

from typing import Callable, Tuple, List, Sequence, Dict
from NumericalMethodSolver import NumericalMethodSolver
from IntermediateValueTheorem import sign
import random

class SecantSolver(NumericalMethodSolver):
    @staticmethod
    def calculate_next_value(func : Callable, x_0 : float, x_1 : float, check_for_root : bool = False) -> float:
        f_x_0, f_x_1 = func(x_0), func(x_1)
        x_new = x_1 - f_x_1 * (x_1 - x_0) / (f_x_1 - f_x_0)

        # Provision for Regula-Falsi Method
        if check_for_root:
            f_x_new = func(x_new)
            return x_new, f_x_new, sign(f_x_new * f_x_1) == -1
        else:
            return x_new
        
    def set_roots(self, init_approx : Sequence[float] = None) -> None:
        if init_approx == None:
            assert not self.approximations == None, "Approximations need to be set before calculating roots by calling self.set_approximations(*args, **kwargs)"
            self.roots = []
            for approximation in self.approximations.keys():
                x_0 = approximation
                x_1 = approximation + 1e-01 * random.random()
                f_x_k = self.func(x_1)
                while(abs(x_1 - x_0) >= self.epsilon and abs(f_x_k) >= self.epsilon):
                    x_new = self.calculate_next_value(self.func, x_0, x_1)
                    x_0 = x_1
                    x_1 = x_new
                    f_x_k = self.func(x_new)
                    
                self.roots.append(x_1)
        else:
            self.roots = []
            if isinstance(init_approx, float) or isinstance(init_approx, int):
                init_approx = [init_approx]
            for i in range(len(init_approx)):
                x_0 = init_approx[i]
                x_1 = init_approx[i] + 1e-01 * random.random()
                f_x_k = self.func(x_1)
                while(abs(x_1 - x_0) >= self.epsilon and abs(f_x_k) >= self.epsilon):
                    x_new = self.calculate_next_value(self.func, x_0, x_1)
                    x_0 = x_1
                    x_1 = x_new
                    f_x_k = self.func(x_new)
                    
                self.roots.append(x_1)
                
