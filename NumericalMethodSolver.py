# A generic Solver Interface that needs to be implemented by any class which aims at solving Non-linear Algaebric Equations for helper functions
# In order to extend the interface, the abstract method set_roots(self, *args, **kwargs) needs to be implemented

from typing import Dict, List, Callable, Sequence
from Bisection import bisect_interval

class NumericalMethodSolver():
    def __init__(self, func : Callable, intervals : Sequence = None, epsilon : float = 1e-5) -> None:
        self.func = func
        self.intervals = intervals
        self.epsilon = epsilon
        self.approximations = None
        self.roots = None
    
    def get_approximations(self) -> Dict:
        assert not self.approximations == None, "Approximations need to be set by calling self.set_approximations(*args, **kwargs)"
        return self.approximations

    def get_roots(self) -> List:
        assert not self.roots == None, "Roots need to be calculated by calling self.set_roots(*args, **kwargs)"
        return self.roots

    def set_intervals(intervals:  Sequence) ->None:
        self.intervals = intervals
    
    def set_approximations(self, iterations : int = None, epsilon : float = None, digits : int = 5) -> None:
        assert not self.intervals == None, "Intervals need to be set for estimating approximations using self.set_intervals(*args, **kwargs)"
        self.approximations = bisect_interval(self.func, self.intervals, iterations, epsilon, digits)
        
    def set_roots(self, *args, **kwargs) -> None:
        pass

    def print_roots(self, *args, **kwargs) -> None:
        assert not self.roots == None, "Roots need to be calculated by calling self.set_roots(*args, **kwargs)"
        print(self.roots) 
