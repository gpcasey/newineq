__author__ = 'Greg'
import InequalityClasses as ineq
import matplotlib as plot



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
A_0 = 10.0
B_RICH_0 = .001 / RICH_FRAC
B_POOR_0 = 0.0
EDU_POOR_0 = 0.0
EDU_RICH_0 = 0.0
assert B_RICH_0 >= EDU_RICH_0
assert B_POOR_0 >= EDU_POOR_0


####################################Run Model

#Initialize
first_poor = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_POOR_0, educ=EDU_POOR_0, type="p")
first_rich = ineq.Individual(beta = BETA, theta_bar=THETA_BAR, bequest=B_RICH_0, educ=EDU_RICH_0, type="r")

first_school = ineq.Education(gamma=GAMMA, tech=EDU_TECH)

K_0 = RICH_FRAC * (B_RICH_0 - EDU_RICH_0) + (1 - RICH_FRAC) * (B_POOR_0 - EDU_POOR_0)
print "K_0: ", K_0
H_0 = RICH_FRAC * first_school.hc(EDU_RICH_0) + (1 - RICH_FRAC) * first_school.hc(EDU_POOR_0)
print "H_0: ", H_0

first_firm = ineq.CESProd(sigma=SIGMA, epsilon=EPSILON, K=K_0, H=H_0, A=A_0, delta=DELTA)

current_economy = ineq.Economy(indiv_r=first_rich, indiv_p=first_poor, school=first_school,
                               firm=first_firm, fraction_r=RICH_FRAC, imperfections=IMPERFECTIONS)



#Assertions for getting through all stages
#print current_economy.no_regimeI_ss(krange=KRANGE, kstep=KSTEP)


#Run Economy
inequality_list = []
capital_share_list = []
k_list = []
K_list = []
H_list = []
constrained_list = []
zero_list = []
count = 0
for t in range(MAXT):
    inequality_list.append(current_economy.get_inequality())
    capital_share_list.append(current_economy.get_Kshare())
    k_list.append(current_economy.get_littlek())
    K_list.append(current_economy._firm._K)
    constrained_list.append(current_economy.poor_constrained())
    zero_list.append(current_economy.poor_zero())
    new_econ = current_economy.start_next(krange=KRANGE, kstep=KSTEP)
    current_economy = new_econ
    count += 1
    print count



########################################Print Output
print "Zero", zero_list
print "CONSTRAINT", constrained_list
print "k", k_list
print "K", K_list
print "KSHARE", capital_share_list
print "INEQUALITY", inequality_list
