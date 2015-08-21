__author__ = 'Greg'
import InequalityClasses as ineq

test_school = ineq.Education(.75, 2)

def test_education():
    if not test_school.hc(0) == 1:
        print "EDUC1", test_school.hc(0)
    if not round(test_school.hc(4), 2) == 6.66:
        print "EDUC2", test_school.hc(4)
    if not round(test_school.hprime(3), 2) == 1.14:
        print "EDU3", round(test_school.hprime(3), 2)
    if not round(test_school.hprime_inverse(3), 2) == .06:
        print "EDU4", round(test_school.hprime_inverse(3), 2)

test_education()
