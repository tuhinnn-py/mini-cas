# TODO: Use the Muller's Method for finding out the real roots of a function f.
# Note: Muller's Method has a superlinear order of convergence (1.8)
# Note: Muller's Method is used in cases where f' or f'' doesnot exist.

# Disclaimer:
# Muller's Method is sensitive to initial approximations unlike Newton-Rhapson's Method or Secant Method, since while solving quadratic equations
# the value of the discriminant might end up negative leading to a math domain error, or even worse the value of the determinant might end up to be 0
# leading to a zero-division error. To account for such errors, restart iterations with the acquired roots as the initial approximations to the roots.

# Muller's Method can be proved in the following way:
#   When the function meets the X-axis, the function can be approximated with a quadratic curve(parabola) in an infinitesimally small neighborhood (a2.x^2 + a1.x + a0)
#   A quadratic curve can be written as f(x) = a2.x^2 + a1.x + a0
#   Given, x_0, x_1 and x_2 are the initial approximations to the root,
#   For the purpose of the proof, f(x) can also be written as a2.(x - x_2)^2 + a1.(x - x_2) + a0, where x_2 is the third initial approximation to the root
#   f(x) = a2.(x - x_2)^2 + a1.(x - x_2) + a0 = 0           -(0)
#   f(x_2) = a0                                             -(1)
#   Using (1),
#   f(x_1) - f(x_2) = a2.(x_1 - x_2)^2 + a1.(x_1 - x_2)     -(2)
#   f(x_0) - f(x_2) = a2.(x_0 - x_2)^2 + a1.(x_0 - x_2)     -(3)
#   Solving the system of linear equations (2) and (3) using Cramer's rule,
#   a2 = ((f(x_1) - f(x_2)).(x_0 - x_2) - (f(x_0) - f(x_2)).(x_1 - x_2)) / D
#   a1 = ((f(x_0) - f(x_2)).(x_1 - x_2)^2 - (f(x_1) - f(x_2)).(x_0 - x_2)^2) / D
#   where, D = det|(x_1 - x_2)^2    (x_1 - x_2)| = (x_2 - x_0)(x_2 - x_1)(x_1 - x_0)
#                 |(x_0 - x_2)^2    (x_0 - x_2)|
#   Using Shreedharacharya's Method for (0) and placing the values of the coefficients,
#   x_new = x_2 + (-a1 +- (a1^2 - 4.a2.a0)^(1/2))/(2.a2)
#   Discard the solution farthest away from x_2

# Note: Stop iteration in one of two situations:
#           |f(x)| < epsilon
#           |x(k + 1) - x(k)| < epsilon

# Note: Set manual approximations while calculating roots of f using init_approx
# Note: Instead of using manual approximations, set approximations of f using Bisection Method through self.set_approximations(*args, **kwargs)

from typing import Callable, Tuple, List, Sequence, Dict
from NumericalMethodSolver import NumericalMethodSolver
from IntermediateValueTheorem import sign
import math
import random

class MullerSolver(NumericalMethodSolver):
    @staticmethod
    def get_det(arr : List[List[float]]) -> float:
        return arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0]

    @staticmethod
    def calculate_roots_sd(a : float, b : float, c : float) -> float:
        # Since, we discard the root farthest away from the current iteration value, we need not solve the entire equation
        # Rather, we choose the denominator of the rationalized Shreedharacharya's equation depending on the sign of b
        # x* = x - (1/2).(b -+ (b*b - 4*a*c)^(1/2))/a
        # Upon rationalizing the equation, x* = x - 2*c/(b +- (b*b - 4*a*c)^(1/2))

        try:
            d = math.sqrt(b*b - 4*a*c)
        except:
            print("Math Domain Error. Manually setting discriminant to zero...")
            d = 1e-01
            
        if sign(b) >= 0:
            return - 2 * c / (b + d)
        else:
            return - 2 * c / (b - d)
    
    def set_roots(self, init_approx : float = None) -> None:
        det_is_zero = False
        if init_approx == None:
            assert not self.approximations == None, "Approximations need to be set before calculating roots by calling self.set_approximations(*args, **kwargs)"
            self.roots = []
            for approximation in self.approximations.keys():
                x_0 = approximation
                x_1 = x_0 + 1e-01 * random.random()
                x_2 = x_0 + 1e-01 * random.random()
                
                a_0 = self.func(x_2)
                f_x_0 = self.func(x_0)
                f_x_1 = self.func(x_1)
                
                while(abs(a_0) >= self.epsilon):
                    det = (x_2 - x_0) * (x_2 - x_1) * (x_1 - x_0)
                    
                    # Create a matrix for ease of determinant calculation
                    iter_matrix = [[f_x_0 - a_0, f_x_1 - a_0], [(x_0 - x_2)**2, (x_1 - x_2)**2]]
                    if det == 0:
                        det_is_zero = True
                        break
                    a_1 = self.get_det(iter_matrix) / det

                    iter_matrix = [[f_x_1 - a_0, f_x_0 - a_0], [x_1 - x_2, x_0 - x_2]]
                    a_2 = self.get_det(iter_matrix) / det

                    x_new = x_2 + self.calculate_roots_sd(a_2, a_1, a_0)
                    a_0 = self.func(x_new)

                    # Swap variables for the next iteration
                    x_0 = x_1
                    x_1 = x_2
                    x_2 = x_new
                    
                self.roots.append(x_2)
        else:
            self.roots = []
            if isinstance(init_approx, float) or isinstance(init_approx, int):
                init_approx = [init_approx]
            for i in range(len(init_approx)):
                x_0 = init_approx[i]
                x_1 = x_0 + 1e-01 * random.random()
                x_2 = x_0 + 1e-01 * random.random()
                
                a_0 = self.func(x_2)
                f_x_0 = self.func(x_0)
                f_x_1 = self.func(x_1)
                
                while(abs(a_0) >= self.epsilon):
                    det = (x_2 - x_0) * (x_2 - x_1) * (x_1 - x_0)
                    
                    # Create a matrix for ease of determinant calculation
                    iter_matrix = [[f_x_0 - a_0, f_x_1 - a_0], [(x_0 - x_2)**2, (x_1 - x_2)**2]]
                    if det == 0:
                        det_is_zero = True
                        break
                    a_1 = self.get_det(iter_matrix) / det

                    iter_matrix = [[f_x_1 - a_0, f_x_0 - a_0], [x_1 - x_2, x_0 - x_2]]
                    a_2 = self.get_det(iter_matrix) / det

                    x_new = x_2 + self.calculate_roots_sd(a_2, a_1, a_0)
                    a_0 = self.func(x_new)

                    # Swap variables for the next iteration
                    x_0 = x_1
                    x_1 = x_2
                    x_2 = x_new
                   
                self.roots.append(x_2)
                
        if det_is_zero:
            print("Zero Division Error. Restarting iterations ...")
            init_approx = self.roots
            self.roots = []
            self.set_roots(init_approx)
