# TODO: Demonstrate Intermediate Value Theorem
# Theorem: If a function f satisfies the condition f(a)f(b) < 0 where a < b, then the function f has (atleast) one or a odd number of real roots in the closed
# interval [a, b]
from typing import Tuple, Callable
sign = lambda x: (x > 0) - (x < 0)

# Check whether either interval is a zero of the function and return None, root
# Check whether the interval satisfies the Intermediate Value Theorem 

def has_root(func : Callable, interval : Tuple) -> bool:
    assert len(interval) == 2, "Length of tuple should be equal to 2"
    if func(interval[0]) == 0:
        return None, interval[0]
    if func(interval[1]) == 0:
        return None, interval[1]
    return sign(func(interval[0]) * func(interval[1])) == -1, None
