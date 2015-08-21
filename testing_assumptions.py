__author__ = 'Greg'
import InequalityClasses as ineq
"""
This file tests the methods used to find a steady states in the economy.
"""



MAXT = 10
KRANGE = 10
KSTEP = 200.0

#Parameters
SIGMA = 1.0
EPSILON = .4
DELTA = 1.0
EDU_TECH = 5
GAMMA = .5
BETA = .25
THETA_BAR = .1
RICH_FRAC = .01
IMPERFECTIONS = True

#Starting Values
A_0 = 10.0
B_RICH_0 = .01 / RICH_FRAC
B_POOR_0 = 0.0
EDU_POOR_0 = 0.0
EDU_RICH_0 = 0.0
assert B_RICH_0 >= EDU_RICH_0
assert B_POOR_0 >= EDU_POOR_0



def check_inverse_littlef():
    new_firm = ineq.CESProd(K=.03, H=1, A=10.0, epsilon=.4, sigma=2)
    if not round(new_firm.inverse_wages(), 2) == 12.25:
        print "IF1", round(new_firm.inverse_wages(), 2)

def check_inverse_wages():
    new_firm = ineq.CESProd(K=.03, H=1, A=10.0, epsilon=.4, sigma=1)
    if not round(new_firm.inverse_wages(.03), 5) == .00056:
        "IW1", round(new_firm.inverse_wages(.03), 5)
    new_firm_B = ineq.CESProd(K=new_firm.inverse_wages(.03), H=1, A=10.0, epsilon=.4, sigma=1)
    if not abs(new_firm_B.factor_payments()[1] - .03) < .001:
        print "IW2", new_firm_B.factor_payments()[1]
    new_firm_C = ineq.CESProd(K=0.02098525515420107, H=1, A=10.0, epsilon=.4, sigma=1)
    if not new_firm_C.factor_payments()[1] > new_firm_B.factor_payments()[1]:
        print "IW3", new_firm_C.factor_payments()[1], new_firm_B.factor_payments()[1]
    print new_firm_C.factor_payments()[1], new_firm_B.factor_payments()[1]

def check_khat():
    first_poor = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_POOR_0, educ=EDU_POOR_0, type="p")
    first_rich = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_RICH_0, educ=EDU_RICH_0, type="r")
    first_school = ineq.Education(gamma=GAMMA, tech=EDU_TECH)
    first_firm = ineq.CESProd(K=.03, H=1, A=10.0, epsilon=.4, sigma=1)
    current_economy = ineq.Economy(indiv_r=first_rich, indiv_p=first_poor, school=first_school,
                               firm=first_firm, fraction_r=RICH_FRAC, imperfections=IMPERFECTIONS)

    if not round(current_economy.find_khat(), 5) == .00056:
        print "KHAT1", round(current_economy.find_khat(), 5)


def check_nextk_regimeI():
    pass

check_inverse_wages()
check_khat()