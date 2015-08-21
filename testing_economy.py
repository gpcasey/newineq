__author__ = 'Greg'
import InequalityClasses as ineq



test_poor_c = ineq.Individual(beta = .25, theta_bar = .5, bequest = 3,educ = 3, type = "p")
test_rich = ineq.Individual(beta = .25, theta_bar = .5, bequest = 6, educ = 4, type = "r")
test_school = ineq.Education(gamma=.9, tech=2)
test_firm = ineq.CESProd(sigma = 3, epsilon = .3, A = 2, K = 4, H = 2, delta = 1)
test_firm_B = ineq.CESProd(sigma = 4, epsilon = .3, A = 2, K = 6, H = 2, delta = 1)
test_economy_A = ineq.Economy(indiv_p = test_poor_c, indiv_r=test_rich, firm=test_firm, school=test_school, fraction_r=.1)


def older_check():
    """
    checks input to economy functions, which are already checked.
    """
    if not round(test_firm.output(), 2) == 5.1:
       print "OLD1", round(test_firm.output(), 2)
    if not test_firm_B.output() == test_firm_B._H * test_firm_B.littlef():
        print "OLD2", test_firm_B.output(), test_firm_B._H * test_firm_B.littlef()
    if not round(test_firm_B.output(), 2) == 6.17:
        print "OLD2b", test_firm_B.output()
    if not (test_firm_B.factor_payments()[0] / test_firm_B.factor_payments()[1]
                - ((test_firm_B.get_epsilon() / (1 - test_firm_B.get_epsilon()))
                        * (test_firm_B.get_littlek() ** ( -1 / test_firm_B.get_sigma())))
                            < .001):
        print "OLD3", test_firm_B.factor_payments()
    rho_a, wage_a = test_firm.factor_payments()
    if not (round(rho_a, 3), round(wage_a, 3)) == (round(.516, 3), round(1.518, 3)):
        print "OLD4", (round(rho_a, 2), round(wage_a, 2))
    hc_a = test_school.hc(test_poor_c.get_edu())
    if not round(hc_a, 2) == round(6.38, 2):
        print "OLD5", hc_a


def checking_combinations():
    """
    Checks definitions other than finding k.
    """

    if not round(test_economy_A.get_income(test_poor_c), 2)  == round(1.518 * 6.376, 2):
         print "comb1", round(test_economy_A.get_income(test_poor_c), 2), round(1.518 * 6.376, 2)
    if not test_economy_A.get_income(test_rich)  - (2 * .5197 + 7.96 * 1.528) < .01:
         print "comb2", test_economy_A.get_income(test_rich), (2 * .5197 + 7.96 * 1.518)
    if not round(test_economy_A.get_inequality(), 2) == round((2 * .5197 + 7.96 * 1.518) / (1.518 * 6.376), 2):
         print "comb4", test_economy_A.get_inequality(), round(2 * .5197 + 7.96 * 1.518)
    if not round(test_economy_A.optimal_e(3), 2) == round(66511552.68, 2):
         print "comb5", test_economy_A.optimal_e(3), round(66511552.68, 2)

def checking_bdiv():
    new_firm = ineq.CESProd(A=1, sigma=1, K=1, H = 1, epsilon=.6, delta=1)
    new_school = ineq.Education(gamma=.5, tech=1)
    poor1 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=0, type="p")
    poor2 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=10, type="p")
    rich = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=100, type="r")

    new_economy = ineq.Economy(poor1, rich, new_school, new_firm, .5, imperfections=True)

    # Tests with poor leaving no bequests
    if not new_economy.optimal_e(3) == 1:
        print "BDIV1", new_economy.optimal_e(3)
    if not new_economy.bequest_division(poor1, 3) == (0, 0, 1, 0):
        print "BDIV2", new_economy.bequest_division(poor1, 3)
    if not ([round(new_economy.bequest_division(rich, 3)[0], 2), round(new_economy.bequest_division(rich, 3)[1], 2)]
            == [14.35, 13.35]):
        print "BDIV2b", new_economy.bequest_division(rich, 3)
    if not ([round(new_economy.bequest_division(rich, 3)[2], 2), round(new_economy.bequest_division(rich, 3)[3], 2)]
            == [2, 1]):
        print "BDIV2c", new_economy.bequest_division(rich, 3)

    # poor leave but are constrained.
    new_economy_B = ineq.Economy(poor2, rich, new_school, new_firm, .5, imperfections=True)
    if not new_economy_B.bequest_division(poor2, 3)[0] == new_economy_B.bequest_division(poor2, 3)[3]:
        print "BDIV3a", new_economy_B.bequest_division(poor2, 3)[0], new_economy_B.bequest_division(poor2, 3)[3]
    if not new_economy_B.bequest_division(poor2, 3)[1] == 0:
        print "BDIV3b", new_economy_B.bequest_division(poor2, 3)[1]

    #Poor are constrained but there are no credit market imperfections.
    new_economy_C = ineq.Economy(poor2, rich, new_school, new_firm, .5, imperfections=False)
    if not new_economy_C.bequest_division(poor2, 3)[3] == new_economy_B.optimal_e(3):
        print "BDIV4", new_economy_C.bequest_division(poor2, 3)[3], new_economy_B.optimal_e(3)
    if not new_economy_C.bequest_division(poor2, 3)[1] - (.85 - new_economy_B.optimal_e(3) < .001):
        print "BDIV5", new_economy_C.bequest_division(poor2, 3)[1], .85 - new_economy_B.optimal_e(3)


def check_calc_littlek():
    new_firm = ineq.CESProd(A=1, sigma=1, K=1, H = 1, epsilon=.6, delta=1)
    new_school = ineq.Education(gamma=.5, tech=1)
    poor1 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=0, type="p")
    poor2 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=5, type="p")
    rich = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=100, type="r")

    #Poor dont leave bequests
    new_economy = ineq.Economy(poor1, rich, new_school, new_firm, .5, imperfections=True)
    test_k = 7.2
    assert new_economy.bequest_division(poor1, 10)[0] == 0, "poor too rich"
    assert new_economy.optimal_e(test_k) < new_economy.bequest_division(rich, 10)[0], "rich too poor"
    guessed_num = 1.0 * new_economy._fraction_r* (rich._beta * (new_economy.get_income(rich) - rich.calc_theta())
                                            - new_economy.optimal_e(test_k))
    guessed_denom = 1.0 * (new_economy._fraction_r * (new_school.hc(new_economy.optimal_e(test_k)))
                     + (1 - new_economy._fraction_r))
    guessed = guessed_num / guessed_denom
    if not new_economy.calc_littlek(test_k)[2] == guessed:
        print "littlek1", guessed, new_economy.calc_littlek(test_k)[2]

    #Poor do leave bequests. Are constrained.
    new_economy_B = ineq.Economy(poor2, rich, new_school, new_firm, .5, imperfections=True)
    test_k = 4
    assert new_economy_B.optimal_e(test_k) > new_economy_B.bequest_division(poor2, 10)[0], "poor too rich"
    assert new_economy_B.bequest_division(poor2, 10)[0] > 0, "poor too poor"
    assert new_economy_B.optimal_e(test_k) < new_economy_B.bequest_division(rich, 10)[0], "rich too poor"
    guessed_num = 1.0 * (rich._beta * ((new_economy_B.get_income(rich) * new_economy_B._fraction_r
                                        + new_economy_B.get_income(poor2) * (1 -new_economy_B._fraction_r ))
                                       - rich.calc_theta())
                    - (new_economy_B._fraction_r * new_economy_B.optimal_e(test_k))
                         - ((1 - new_economy_B._fraction_r) * new_economy_B.bequest_division(poor2, test_k)[0]))
    guessed_denom = 1.0 * (new_economy_B._fraction_r * new_school.hc(new_economy_B.optimal_e(test_k))
                     + (1 - new_economy_B._fraction_r) * new_school.hc(new_economy_B.bequest_division(poor2, test_k)[0]))
    #NOTE: This doesnt work with total output becuase there is no assumption that "new_economy" is in equilibrium.
    #       It might make sense to go back and do a functional check for this afterwards.
    guessed = guessed_num / guessed_denom
    if not new_economy_B.calc_littlek(test_k)[2] == guessed:
        print "littlek2", guessed, new_economy_B.calc_littlek(test_k)[2]


    # Poor have optimal level.
    test_k = .000001
    assert new_economy_B.optimal_e(test_k) < new_economy_B.bequest_division(poor2, 10)[0], "poor too poor"
    assert new_economy_B.optimal_e(test_k) < new_economy_B.bequest_division(rich, 10)[0], "rich too poor"
    guessed_num = 1.0 * (rich._beta * ((new_economy_B.get_income(rich) * new_economy_B._fraction_r
                                        + new_economy_B.get_income(poor2) * (1 -new_economy_B._fraction_r ))
                                       - rich.calc_theta())
                    - new_economy_B.optimal_e(test_k))
    guessed_denom = 1.0 * new_school.hc(new_economy_B.optimal_e(test_k))
    #NOTE: This doesnt work with total output becuase there is no assumption that "new_economy" is in equilibrium.
    #       It might make sense to go back and do a functional check for this afterwards.
    guessed = guessed_num / guessed_denom
    if not new_economy_B.calc_littlek(test_k)[2] == guessed:
        print "littlek3", guessed, new_economy_B.calc_littlek(test_k)[2]


    # Poor are contrained, but no credit market imperfections
    new_economy_B = ineq.Economy(poor2, rich, new_school, new_firm, .5, imperfections=False)
    test_k = 4
    assert new_economy_B.optimal_e(test_k) > new_economy_B.bequest_division(poor2, 10)[0], "poor too rich"
    assert new_economy_B.bequest_division(poor2, 10)[0] > 0, "poor too poor"
    assert new_economy_B.optimal_e(test_k) < new_economy_B.bequest_division(rich, 10)[0], "rich too poor"
    guessed_num = 1.0 * (rich._beta * (new_economy_B.get_income(rich) * new_economy_B._fraction_r
                    + new_economy_B.get_income(poor2) * (1 -new_economy_B._fraction_r ) - rich.calc_theta())
                         - new_economy_B.optimal_e(test_k))
    guessed_num_B =  ((1 - new_economy_B._fraction_r) * new_economy_B.bequest_division(poor2, test_k)[1]
        + new_economy_B._fraction_r * new_economy_B.bequest_division(rich, test_k)[1])
    guessed_denom = 1.0 * new_school.hc(new_economy_B.optimal_e(test_k))

    #NOTE: This doesnt work with total output becuase there is no assumption that "new_economy" is in equilibrium.
    #       It might make sense to go back and do a functional check for this afterwards.
    guessed = guessed_num / guessed_denom
    if not new_economy_B.calc_littlek(test_k)[2] == guessed:
        print "littlek4", guessed, new_economy_B.calc_littlek(test_k)[2]
        print guessed_num_B,  new_economy_B.bequest_division(poor2, test_k)[1],  new_economy_B.bequest_division(rich, test_k)[1]
        print "K", guessed_num, new_economy_B.calc_littlek(test_k)[0]
        print "H", guessed_denom, new_economy_B.calc_littlek(test_k)[1]


def check_pfe():
    """
    Checks finding of optimal k. Shows minimization.
    """
    new_firm = ineq.CESProd(A=1, sigma=1, K=1, H = 1, epsilon=.6, delta=1)
    new_school = ineq.Education(gamma=.5, tech=1)
    poor1 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=0, type="p")
    poor2 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=5, type="p")
    poor3 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=100, type="p")
    rich = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=100, type="r")

    # Poor dont leave bequests
    new_economy = ineq.Economy(poor1, rich, new_school, new_firm, .5, imperfections=True)
    newk = new_economy.find_nextk(10, 1000.0)[0]
    if not abs(newk - new_economy.calc_littlek(newk)[2]) < .001:
        print "pfe1", new_economy.find_nextk(10, 1000.0), new_economy.calc_littlek(newk)[2]
    assert new_economy.bequest_division(poor1, newk)[0] == 0

    #Poor do leave bequests. Are constrained.
    new_economy_B = ineq.Economy(poor2, rich, new_school, new_firm, .5, imperfections=True)
    newk = new_economy_B.find_nextk(10, 1000.0)[0]
    if not abs(newk - new_economy_B.calc_littlek(newk)[2]) < .001:
        print "pfe2", new_economy_B.find_nextk(10, 1000.0), new_economy_B.calc_littlek(newk)[2]
    assert new_economy_B.bequest_division(poor2, newk)[0] < new_economy_B.optimal_e(newk)

    #Poor do leave bequests. Not constrained.
    new_economy_C = ineq.Economy(poor3, rich, new_school, new_firm, .5, imperfections=True)
    newk = new_economy_C.find_nextk(5, 2000.0)[0]
    if not abs(newk - new_economy_C.calc_littlek(newk)[2]) < .001:
        print "pfe3", new_economy_C.find_nextk(10, 1000.0), new_economy_C.calc_littlek(newk)[2]
    assert new_economy_C.bequest_division(poor3, newk)[0] > new_economy_C.optimal_e(newk)

    # No Credit market imperfections
    new_economy_D = ineq.Economy(poor2, rich, new_school, new_firm, .5, imperfections=False)
    newk = new_economy_D.find_nextk(10, 1000.0)[0]
    if not abs(newk - new_economy_D.calc_littlek(newk)[2]) < .001:
        print "pfe2", new_economy_D.find_nextk(10, 1000.0), new_economy_D.calc_littlek(newk)[2]
    assert new_economy_D.bequest_division(poor2, newk)[1] < 0


def checking_startnext():
    new_firm = ineq.CESProd(A=1, sigma=1, K=1, H = 1, epsilon=.6, delta=1)
    new_school = ineq.Education(gamma=.5, tech=1)
    poor1 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=0, type="p")
    poor2 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=5, type="p")
    poor3 = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=100, type="p")
    rich = ineq.Individual(beta = .25, theta_bar=1, educ=0, bequest=100, type="r")

    # Poor dont leave bequests
    new_economy = ineq.Economy(poor1, rich, new_school, new_firm, .5, imperfections=True)
    newk = new_economy.find_nextk(10, 1000.0)[0]
    second_economy = new_economy.start_next(10, 1000.0)
    if not round(second_economy.get_littlek(),2) == round(newk, 2):
        print "NEXT1", second_economy.get_littlek(), newk
    # print new_economy.get_littlek()
    # print second_economy.get_littlek()
    # current_econ = second_economy
    # for dummy_i in range(10):
    #     next_econ =  current_econ.start_next(10, 1000.0)
    #     print next_econ.get_littlek()

    #Poor do leave bequests. Are constrained.
    new_economy_B = ineq.Economy(poor2, rich, new_school, new_firm, .5, imperfections=True)
    newk = new_economy_B.find_nextk(10, 1000.0)[0]
    second_economy = new_economy_B.start_next(10, 1000.0)
    if not round(second_economy.get_littlek(),2) == round(newk, 2):
        print "NEXT2", second_economy.get_littlek(), newk

    #Poor do leave bequests. Not constrained.
    new_economy_C = ineq.Economy(poor3, rich, new_school, new_firm, .5, imperfections=True)
    newk = new_economy_C.find_nextk(5, 2000.0)[0]
    second_economy = new_economy_C.start_next(5, 2000.0)
    if not round(second_economy.get_littlek(),2) == round(newk, 2):
        print "NEXT3", second_economy.get_littlek(), newk

    # No Credit market imperfections
    new_economy_D = ineq.Economy(poor2, rich, new_school, new_firm, .5, imperfections=False)
    newk = new_economy_D.find_nextk(10, 1000.0)[0]
    second_economy = new_economy_D.start_next(10, 1000.0)
    if not round(second_economy.get_littlek(),2) == round(newk, 2):
        print "NEXT4", second_economy.get_littlek(), newk


older_check()
checking_combinations()
checking_bdiv()
check_calc_littlek()
check_pfe()
print "**********************Start Next********************"
checking_startnext()