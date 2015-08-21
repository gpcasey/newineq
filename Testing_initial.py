__author__ = 'Greg'
import InequalityClasses as ineq
import matplotlib as plot



#####################################Define Constants

# Technical Constants
MAXT = 10
KRANGE = 50
KSTEP = 1000.0

#Parameters
SIGMA = 1.0
EPSILON = .4
DELTA = 1.0
EDU_TECH = 5
GAMMA = .25
BETA = .25
THETA_BAR = .25
RICH_FRAC = .01
IMPERFECTIONS = True

#Starting Values
A_0 = 2.0
B_RICH_0 = .001 / .01
B_POOR_0 = 0.0
EDU_POOR_0 = 0.0
EDU_RICH_0 = 0.0
K_0 = RICH_FRAC * (B_RICH_0 - EDU_RICH_0) + (1 - RICH_FRAC) * (B_POOR_0 - EDU_POOR_0)
print "K_0: ", K_0

####################################Run Model

#Initialize
first_poor = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_POOR_0, educ=EDU_POOR_0, type="p")
first_rich = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_RICH_0, educ=EDU_RICH_0, type="r")

first_school = ineq.Education(gamma=GAMMA, tech=EDU_TECH)
H_0 = RICH_FRAC * first_school.hc(EDU_RICH_0) + (1 - RICH_FRAC) * first_school.hc(EDU_POOR_0)
print "H_0: ", H_0

first_firm = ineq.CESProd(sigma=SIGMA, epsilon=EPSILON, K=K_0, H=H_0, A=A_0, delta=DELTA)

current_economy = ineq.Economy(indiv_r=first_rich, indiv_p=first_poor, school=first_school,
                               firm=first_firm, fraction_r=RICH_FRAC, imperfections=IMPERFECTIONS)


print current_economy.bequest_division(first_rich, 0.07894074743455817)