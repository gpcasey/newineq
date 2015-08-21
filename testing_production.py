__author__ = 'Greg'
import InequalityClasses as ineq

test_firm_ces = ineq.CESProd(2, .25, 64, 16)
test_firm_cd = ineq.CESProd(1, .25, 64, 16)
test_firm_cesB = ineq.CESProd(2, .25, 64, 16, A = 2) # A added later

def check_output():
    """
    Check output calculations
    """
    #Big F
    if not  test_firm_ces.output() == 25:
        print "OUT1", test_firm_ces.output()
    if not abs(test_firm_cd.output() - 22.63) < .01:
        print "OUT2", test_firm_cd.output()

    #little F (and equality)
    if not  test_firm_ces.littlef() == 25 / 16.0:
        print "OUT3", test_firm_ces.littlef(), 25 / 16
    if not abs(test_firm_cd.littlef() - 22.63 / 16) < .001:
        print "OUT4", test_firm_cd.output()

    # A test (Added later)
    test_firm_cesB = ineq.CESProd(2, .25, 64, 16, A = 2)
    if not  test_firm_cesB.output() == test_firm_ces.output() * 2:
        print "OUT5", test_firm_ces.output()
    if not test_firm_cesB.littlef() == test_firm_ces.littlef() * 2:
        print "OUT6", test_firm_cesB.littlef()


def check_factor_payments():
    """
    check factor payments calculations.
    NOTE:  delta not tested.
    """

    # Check Calculation.
    bigR, wages = test_firm_ces.factor_payments()
    if not max([abs(bigR - .15625 ), abs(wages - .9375)]) < .001:
        print "FAC1", (bigR, wages)

    bigR, wages = test_firm_cd.factor_payments()
    if not max([abs(bigR - .088), abs(wages - 1.061)]) < .001:
        print "FAC2", (bigR, wages)

    # check A
    bigR, wages = test_firm_ces.factor_payments()
    bigR_B, wages_B = test_firm_cesB.factor_payments()
    if not max([abs(2.0 * bigR - bigR_B), abs(2.0 * wages - wages_B)]) < .00001:
        print "FAC1B", (2.0 * bigR, 2.0 * wages), (bigR_B, wages_B )

    # Check equality for wages with f(k) - f'(k)k
    bigR, wages = test_firm_ces.factor_payments()
    fk = test_firm_ces.littlef()
    littlek = test_firm_ces._K / test_firm_ces._H
    fprimek = test_firm_ces._epsilon * fk ** (1.0 / test_firm_ces._sigma) * littlek ** (-1.0 / test_firm_ces._sigma)
    if not bigR - fprimek == 0.0:
        print "FAC3", bigR - fprimek
    if not wages - (fk - fprimek * littlek) == 0.0:
        print "FAC4", wages - (fk - fprimek * littlek)

    # check limits  (restricted to the case of sigma >= 1)
    limit_test_low = ineq.CESProd(2, .25, .00000000000000000000000000000000000000000000000000000000000000000000000001, .01)
    bigR, wages = limit_test_low.factor_payments()
    if bigR < 100000:
        print "FAC5", bigR


    limit_test_high = ineq.CESProd(2, .25, 10000000000000000000000000000000000000000000000000000000000000000000000000, .01)
    bigR, wages = limit_test_high.factor_payments()
    if bigR != limit_test_low._epsilon ** 2:
        print "FAC6", bigR

    limit_test_high = ineq.CESProd(2, .25, 1, 10000000000000000000000000000000000000000000000000000000000000000000000000)
    bigR, wages = limit_test_high.factor_payments()
    if wages != (1 - limit_test_low._epsilon) ** 2:
        print "FAC7", wages

    limit_test_high = ineq.CESProd(2, .25, 1, .00000000000000000000000000000000000000000000000000000000000000000000000001)
    bigR, wages = limit_test_high.factor_payments()
    if wages < 100000000000:
        print "FAC8", wages

def check_capital_share():
    if not test_firm_ces.capital_share() == 1.6 * test_firm_ces._epsilon:
        print "CSHARE1", test_firm_ces.capital_share()
    if not abs(test_firm_cd.capital_share() - test_firm_cd._epsilon) < .000001:
        print "CSHARE2", test_firm_cd.capital_share()
    if not test_firm_cesB.capital_share() == 1.6 * test_firm_ces._epsilon:
        print "CSHARE3", test_firm_cesB.capital_share()

check_output()
check_factor_payments()
check_capital_share()