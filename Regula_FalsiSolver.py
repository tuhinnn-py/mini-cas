# TODO: Use the Regula-Falsi Method for finding out the real roots of a function f.
# Note: Regula-Falsi Method has a first order of convergence (1)

# Method : c = b - f(b) * (b - a) / (f(b) - f(a)), where a, b are intital approximations to the real root of the function
# Method : Replace a with b and b with c for the next iteration if f(b)f(c) < 1 else Retain a and replace b with c for the next iteration
# Note: Regula-Falsi Method uses two initial approximations to the root

# Regula-Falsi Method can be proved in the following way:
#   When the curve meets the X-axis, the curve can be approximated with a straight line (a1.x + a0) in an infinitesimally small neighborhood
#   Calculate the values of the coefficients(a0, a1) by solving the following system of equations : a1.a + a0 = f(a), a1.b + a0 = f(b): c = -a0 / a1

# Note: Stop iteration in one of two situations:
#           |f(x)| < epsilon
#           |x(k + 1) - x(k)| < epsilon

# Note: Set manual approximations while calculating roots of f using init_approx
# Note: Instead of using manual approximations, set approximations of f using Bisection Method through self.set_approximations(*args, **kwargs)


# Note: Regula-Falsi Method is simply a variation of the Secant Method where we use Intermediate Value Theorem to decide the next interval

from typing import Callable, Tuple, List, Sequence, Dict
from SecantSolver import SecantSolver
import random

class Regula_FalsiSolver(SecantSolver):
    def set_roots(self, init_approx : Sequence[float] = None) -> None:
        if init_approx == None:
            assert not self.approximations == None, "Approximations need to be set before calculating roots by calling self.set_approximations(*args, **kwargs)"
            self.roots = []
            for approximation in self.approximations.keys():
                x_0 = approximation
                x_1 = approximation + 1e-01 * random.random()
                f_x_k = self.func(x_1)
                while(abs(x_1 - x_0) >= self.epsilon and abs(f_x_k) >= self.epsilon):
                    x_new, f_x_k, has_root_ = super(Regula_FalsiSolver, self).calculate_next_value(self.func, x_0, x_1, True)
                    # Adding an if statement simply changes the Secant Solver to a Regula-Falsi Solver
                    if has_root_:
                        x_0 = x_1
                    x_1 = x_new

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
                    x_new, f_x_k, has_root_ = super(Regula_FalsiSolver, self).calculate_next_value(self.func, x_0, x_1, True)
                    # Adding an if statement simply changes the Secant Solver to a Regula-Falsi Solver
                    if has_root_:
                        x_0 = x_1
                    x_1 = x_new

                self.roots.append(x_1)
