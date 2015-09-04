__author__ = 'Greg'
import InequalityClasses as ineq
"""
This file tests the methods used to find a steady states in the economy.
"""



#####################################Define Constants

# Technical Constants
MAXT = 30
KRANGE = 1000
KSTEP = 1000.0
KTOL = .01

#Parameters
SIGMA = 1.5
EPSILON = .7
DELTA = 1.0
EDU_TECH = 1
GAMMA = .5
BETA = .25
THETA_BAR = .01
RICH_FRAC = .001
IMPERFECTIONS = True

#Starting Values
A_0 = 1.0
B_RICH_0 = .001 / RICH_FRAC
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
    new_firm = ineq.CESProd(K=.02, H=1, A=1, epsilon=.5, sigma=1)
    if not round(new_firm.inverse_wages(.03), 5) == round(0.06**2, 5):
        print "IW1", round(new_firm.inverse_wages(.03), 5), round(0.06**2, 5)
    new_firm_B = ineq.CESProd(K=new_firm.inverse_wages(.03), H=new_firm._H, A=new_firm._A, epsilon=new_firm._epsilon, sigma=new_firm.get_sigma())
    if not abs(new_firm_B.factor_payments()[1] - .03) < .001:
        print "IW2", new_firm_B.factor_payments()[1]
    new_firm_C = ineq.CESProd(K=0.02098525515420107, H=new_firm._H, A=new_firm._A, epsilon=new_firm._epsilon, sigma=new_firm.get_sigma())
    if not new_firm_C.factor_payments()[1] > new_firm_B.factor_payments()[1]:
        print "IW3", new_firm_C.factor_payments()[1], new_firm_B.factor_payments()[1]

def check_khat():
    first_poor = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_POOR_0, educ=EDU_POOR_0, type="p")
    first_rich = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_RICH_0, educ=EDU_RICH_0, type="r")
    first_school = ineq.Education(gamma=GAMMA, tech=EDU_TECH)
    first_firm = ineq.CESProd(K=.03, H=1, A=10.0, epsilon=.4, sigma=1)
    current_economy = ineq.Economy(indiv_r=first_rich, indiv_p=first_poor, school=first_school,
                               firm=first_firm, fraction_r=RICH_FRAC, imperfections=IMPERFECTIONS)

    if not round(current_economy.find_khat(), 5) == .00056:
        print "KHAT1", round(current_economy.find_khat(), 5)


def check_regimeI_ss():
    #Parameters
    SIGMA = 1.5
    EPSILON = .7
    DELTA = 1.0
    EDU_TECH = 1
    GAMMA = .5
    BETA = .25
    THETA_BAR = .01
    RICH_FRAC = .001
    IMPERFECTIONS = True

    #Starting Values
    A_0 = 1.0
    B_RICH_0 = .001 / RICH_FRAC
    B_POOR_0 = 0.0
    EDU_POOR_0 = 0.0
    EDU_RICH_0 = 0.0
    assert B_RICH_0 >= EDU_RICH_0
    assert B_POOR_0 >= EDU_POOR_0

    #new parameters
    new_krange = 1          #with these initial conditions the cutoff for regime I should be less than 1.
    new_kstep = 2000.0

    #Initialize economy.
    first_poor = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_POOR_0, educ=EDU_POOR_0, type="p")
    first_rich = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_RICH_0, educ=EDU_RICH_0, type="r")

    first_school = ineq.Education(gamma=GAMMA, tech=EDU_TECH)

    K_0 = RICH_FRAC * (B_RICH_0 - EDU_RICH_0) + (1 - RICH_FRAC) * (B_POOR_0 - EDU_POOR_0)
    H_0 = RICH_FRAC * first_school.hc(EDU_RICH_0) + (1 - RICH_FRAC) * first_school.hc(EDU_POOR_0)

    first_firm = ineq.CESProd(sigma=SIGMA, epsilon=EPSILON, K=K_0, H=H_0, A=A_0, delta=DELTA)

    current_economy = ineq.Economy(indiv_r=first_rich, indiv_p=first_poor, school=first_school,
                                   firm=first_firm, fraction_r=RICH_FRAC, imperfections=IMPERFECTIONS)

    # Find regime I SS
    guess_ss = current_economy.find_ss_regimeI(krange=new_krange, kstep=new_kstep)[0]

    #Start a new economy in the regime I SS to make sure it stays there.
    new_K_0 = guess_ss
    second_rich = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest= guess_ss / RICH_FRAC, educ=EDU_RICH_0, type="r")
    second_firm = ineq.CESProd(sigma=SIGMA, epsilon=EPSILON, K=new_K_0, H=H_0, A=A_0, delta=DELTA)

    second_economy = ineq.Economy(indiv_r=second_rich, indiv_p=first_poor, school=first_school,
                               firm=second_firm, fraction_r=RICH_FRAC, imperfections=IMPERFECTIONS)

    new_economy = second_economy
    for dummy_i in range(30):   #iterate  for 30 periods to make sure still in SS.
        next_economy = new_economy.start_next(krange=new_krange, kstep=new_kstep)
        new_economy = next_economy
    if not abs(new_economy.get_littlek() - guess_ss) < .01:
        print "RIss1", new_economy.get_littlek(), guess_ss


print "******************Inverse Wages******************"
check_inverse_wages()
print "******************khat******************"
check_khat()
print "******************Regime I SS******************"
check_regimeI_ss()
print "**********************No Trivial*****************"
