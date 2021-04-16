import argparse
parser = argparse.ArgumentParser(description = 'Calculate real roots of a polynomial function usinga Numerical Method Solver')
parser.add_argument('solver_', type = str, default = None, help = 'Name of the solver. The system randomly chooses either Secant Solver or Muller Solver by default whenever this program is run')
parser.add_argument('func__', type = function, help = 'Polynomial function whose zeros need to be calculated')
parser.add_argument('der_func__', type = function, default = None, help = 'Derivative of the Polynomial function. 

args = parser.parse_args()
