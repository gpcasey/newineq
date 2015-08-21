__author__ = 'Greg'
import InequalityClasses as ineq

test_indiv_convex = ineq.Individual(.25, 2, 4, 2, "p")
test_indiv_linear = ineq.Individual(.25, 0, 4, 2, "p")



def check_indiv():
    if not test_indiv_convex.income(3, 2, 3) == 12:
        print "IND1", test_indiv.income(3, 2, 3)
    if not test_indiv_linear.decisions(3) == (3 * .75, .25 * 3):
        print "IND2", test_indiv_convex.decisions(3)
    if not test_indiv_convex.decisions(3) == (3, 0):
        print "IND3", test_indiv_convex.decisions(3)
    if not test_indiv_convex.decisions(7) == (6.75, .25):
        print "IND4", test_indiv_convex.decisions(7)

check_indiv()