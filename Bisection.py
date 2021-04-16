# TODO: Bisect a closed interval [a, b] and apply the Intermediate Value Theorem either for a specified number of iterations or untill a particular accuracy
# has already been achieved
# Bisection helps in estimating an accurate initial approximation to a real root of a function, thus making sure that the iterative method does not diverge

from IntermediateValueTheorem import has_root
from typing import Sequence, List, Callable, Tuple
import math

# Calculate the minimum number of iterations required for a given epsilon
def calculate_iterations(interval : Tuple, epsilon : float) -> int:
    return math.ceil(math.log(interval[1] - interval[0], 2) - math.log(epsilon, 2))

# For an accuracy of epsilon e0, the prescribed number of iterations is given by the formula (ln(b - a) - ln(e0)) / ln(2)
# Number of intervals need not necessarily be one, i.e. can be greater than one
# For the approximations, 'digits' number of digits are considered after the decimal place

# Note: For Bisection to work, the half-interval of the given intervals must contain either one or an odd number of real roots
# Half-intervals are defined as follows : (lower-limit, mid-limit) and (mid-limit, upper-limit)

def bisect_interval(func : Callable, intervals : Sequence, iterations : int = None, epsilon : float = None, digits : None = 5) -> List:
    if isinstance(intervals, tuple):
        intervals = [intervals]
    assert intervals[0][0] < intervals[0][1], "Lower limit should be strictly lesser than or equal to Upper limit"
    assert not(iterations == None) or not(epsilon == None), "Either the number of iterations or the degree of accuracy should be specified, Both cannot be None"
    # Domain Constraints
    if not iterations == None:
        assert isinstance(iterations, int) and iterations > 0, "Number of iterations should be an integer and should trivially be greater than zero"
    if not epsilon == None:
        iterations = 0 if iterations == None else iterations
        for interval in intervals:
            iterations = max(iterations, calculate_iterations(interval, epsilon))
    for i in range(iterations):
        new_intervals = []
        while(len(intervals) > 0):
            interval = intervals.pop(0)
            # Calculate the mid-point of the interval
            mid_interval = interval[0] + (interval[1] - interval[0]) / 2
            has_root_, _ = has_root(func, (interval[0], mid_interval))
            if has_root_ or has_root_ == None:
                new_intervals.append((interval[0], mid_interval))

            # If the mid_interval turns out to be a real root, checking the second range is unnecessary
            if not (has_root_ == None and _ == mid_interval):
                has_root_, _ = has_root(func, (mid_interval, interval[1]))
                if has_root_ or has_root_ == None:    
                    new_intervals.append((mid_interval, interval[1]))
        intervals = new_intervals

    # List of initial approximations to the real roots of the function contained in the intervals
    approximations = dict()
    for i in range(len(intervals)):
        approximation = intervals[i][0] + (intervals[i][1] - intervals[i][0]) / 2
        approximations[round(approximation, digits)] = func(approximation)
        
    return approximations

# Perform unit tests with the following examples
# func = lambda x : x**2 - 5*x + 6
# bisect_interval(func, (1, 4), 10) ->  {2.00049: -0.00048804283142089844, 2.99951: -0.00048804283142089844}
# bisect_interval(func, (1, 4), epsilon = 1e-5) ->  {2.0: 9.536752259009518e-07, 3.0: 9.536752259009518e-07}
            
    
        
