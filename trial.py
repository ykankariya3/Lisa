from logaut import ltl2dfa
from pylogics.parsers import parse_ltl

formula = parse_ltl("F(a)")
dfa = ltl2dfa(formula , backend="lydia")


