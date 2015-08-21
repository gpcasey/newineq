__author__ = 'Greg'

class CESProd:
    """
    CES Representative Firm
    """

    def __init__(self, sigma, epsilon, K, H, A = 1, delta = 1):
        """
        Note A is fixed.
        """
        assert sigma >= 1, "weakly increasing capital share."
        if sigma > 1:
            self._nu = 1.0 * sigma / (sigma - 1)
        self._sigma = sigma * 1.0
        self._epsilon = epsilon * 1.0
        self._delta = delta * 1.0
        self._K = K * 1.0
        self._H = H * 1.0
        self._A = A * 1.0

    def get_sigma(self):
        return self._sigma

    def get_epsilon(self):
        return self._epsilon

    def get_delta(self):
        return self._delta

    def get_A(self):
        return self._A

    def get_littlek(self):  #not checked.
        return self._K / self._H

    def output(self):
        """
        Calculate output. CES.
        """
        if self._sigma == 1:
            return self._A * (self._K ** self._epsilon) * (self._H ** (1.0 - self._epsilon))
        else:
            inside = self._epsilon * self._K ** (1.0 / self._nu) + (1 - self._epsilon) * self._H ** (1.0 / self._nu)
            return self._A * (inside ** self._nu)

    def littlef(self):
        """
        f(k). For testing purposes.
        """
        littlek = self._K / self._H
        if self._sigma != 1:
            inside = self._epsilon * littlek ** (1.0 / self._nu) + (1 - self._epsilon)
            return self._A * inside ** self._nu
        else:
            return self._A * littlek ** self._epsilon

    def factor_payments(self):
        """
        Calc return and wages.
        """
        if self._sigma != 1:
            rho = self._A ** (1 - 1.0 / self._sigma) * (self._epsilon * self.output() ** (1.0 / self._sigma)
                             * self._K ** (-1.0 / self._sigma))
            wages = self._A ** (1 - 1.0 / self._sigma) * ((1 - self._epsilon) * self.output() ** (1.0 / self._sigma)
                               * self._H ** (-1.0 / self._sigma))
        else:
            k = self.get_littlek()
            rho = self._A * self._epsilon * k ** (self._epsilon - 1)
            wages = self._A * (1 - self._epsilon) * k ** self._epsilon
        bigR = 1 + rho - self._delta
        return (bigR, wages)


    def inverse_littlef(self, in_y):
        """
        returns k =  f^-1(y).
        """
        if self._sigma != 1:
            in_y = 1.0 * in_y
            num = in_y ** (1 / self._nu) - (1 - self._epsilon)
            denom = self._epsilon
            return (num / denom) ** self._nu
        else:
            return (in_y * 1.0) ** (1 / self._epsilon)



    def inverse_wages(self, in_w):
        """
        returns k = w^-1(wages).
        """
        if self._sigma !=1:
            num = in_w * 1.0
            denom = 1.0 * self._A ** (1 - 1.0 / self._sigma) * (1 - self._epsilon)
            inside = (num / denom) ** self._sigma
            return self.inverse_littlef(inside)
        else:
            num = 1.0 * in_w
            denom = self._A * (1 - self._epsilon)
            inside = num / denom
            return 1.0 * inside ** (1/self._epsilon)

    def capital_share(self):
        """
        returns capital share of income.
        """
        FK, FH = self.factor_payments()
        theta_K = 1.0 * (FK * self._K) / self.output()
        theta_H = 1.0 * (FH * self._H) / self.output()
        assert (theta_K + theta_H)- 1 < .01, "Excess profits"
        return theta_K


class Education:
    """
    Class for the system of education.
    """

    def __init__(self, gamma, tech = 1):
        self._gamma = gamma * 1.0
        self._tech = tech * 1.0

    def copy(self):
        return Education(self._gamma, self._tech)

    def hc(self, educ):
        """
        Calc human capital output.
        """
        return 1 + 1.0 * self._tech * educ ** self._gamma

    def hprime(self, educ):
        """
        return to hc
        """
        return 1.0 * self._gamma * self._tech * educ ** (self._gamma - 1)

    def hprime_inverse(self, hprime):
        """
        inverse of h prime. For calculating e(k).
        """
        frac = 1.0 * hprime / (self._gamma * self._tech)
        return frac ** (1.0 / (self._gamma - 1))


class Individual:
    """
    Class for an individual.
    """
    def __init__(self, beta, theta_bar, bequest, educ, type):
        """
        set-up for individual
        """
        assert type in ("p", "r")
        assert beta > 0 and beta < 1
        self._beta = beta * 1.0
        self._educ = educ * 1.0
        self._bequest = bequest * 1.0
        self._theta_bar = theta_bar * 1.0
        self._type = type

    def get_type(self):
        return self._type

    def get_beta(self):
        return self._beta

    def get_theta_bar(self):
        return self._theta_bar

    def get_bequest(self):
        return self._bequest

    def get_edu(self):
        return self._educ

    def calc_theta(self):
        """
        returns theta
        """
        beta_frac = (1 - self._beta) / self._beta * 1.0
        return self._theta_bar * beta_frac


    def income(self, bigR, wage, hc):
        return 1.0 * wage * hc + bigR * (self._bequest - self._educ)

    def decisions(self, income):
        """
        returns optimal consumption and bequest decisions.
        """
        beta_frac = 1.0 * (1 - self._beta) / self._beta
        theta = self._theta_bar * beta_frac
        if income > theta:
            next_bequest = self._beta * (income - theta)
        else:
            next_bequest = 0
        consumption = income - next_bequest
        return (consumption, next_bequest)


class Economy:
    """
    Class for an economy. Consists of a firms, two indivs, and a school.
    """
    def __init__(self, indiv_p, indiv_r, school, firm, fraction_r, imperfections = False):
        assert isinstance(indiv_p, Individual)
        assert indiv_p.get_type() == "p"
        assert isinstance(indiv_r, Individual)
        assert indiv_r.get_type() == "r"
        assert isinstance(firm, CESProd)
        assert isinstance(school, Education)
        self._firm = firm
        self._school = school
        self._indiv_p = indiv_p
        self._indiv_r = indiv_r
        self._fraction_r = fraction_r
        self._imperfections = imperfections

    def get_Kshare(self):
        return self._firm.capital_share()

    def get_littlek(self):
        return self._firm.get_littlek()

    def poor_constrained(self):
        """
        Returns True if poor are constrained.
        """
        return (self._indiv_p.get_edu() == self._indiv_p.get_bequest())

    def poor_zero(self):
        """
        return True if poor have no bequest
        """
        return self._indiv_p.get_bequest() == 0

    def get_income(self, person):  # could refactor into other functions.
        """
        return income for an individual.
        """
        assert isinstance(person, Individual)
        assert  person in [self._indiv_p, self._indiv_r]
        bigR, wages = self._firm.factor_payments()
        hc = self._school.hc(person.get_edu())
        income = person.income(bigR, wages, hc)
        return income

    def get_inequality(self):
        """
        return inequality.
        """
        poor_income = self.get_income(self._indiv_p)
        rich_income = self.get_income(self._indiv_r)
        inequality = rich_income / poor_income
        assert inequality >= 1
        return inequality

    def optimal_e(self, in_next_k):
        """
        return optimal e(k).
        """
        new_firm = CESProd(sigma=self._firm.get_sigma(), epsilon=self._firm.get_epsilon(), K=in_next_k, H=1,
                               delta = self._firm.get_delta(), A = self._firm.get_A())
        next_bigR, next_wages = new_firm.factor_payments()
        optimal_e = self._school.hprime_inverse(1.0 * next_bigR / next_wages)
        return optimal_e

    def bequest_division(self, indiv, in_next_k):
        """
        Returns bequest and education for the next period.
        """
        assert isinstance(indiv, Individual)
        assert indiv in [self._indiv_p, self._indiv_r]
        optimal_e = self.optimal_e(in_next_k)
        current_income = self.get_income(indiv)
        next_consumption, next_bequest = indiv.decisions(current_income)
        if self._imperfections == False:
            next_edu = optimal_e
        else:
            next_edu = min([optimal_e, next_bequest])
        next_hc = self._school.hc(next_edu)
        next_pc = next_bequest - next_edu  #negative implies borrowing capital.
        return next_bequest, next_pc, next_hc, next_edu


    def calc_littlek(self, in_next_k):
        """
        Helper function for find_nextk().
        Calculates littlek next period given an input level of littlek.
        Does so by computing capital and education investments for each type of person and aggregating.
        """
        next_bigK = 0.0
        next_bigH = 0.0
        for person in [self._indiv_p, self._indiv_r]:
            fraction = (1.0 - 2 * self._fraction_r) * (person == self._indiv_p) + self._fraction_r
            beq, pc, hc, edu = self.bequest_division(person, in_next_k)
            next_bigK += fraction * pc
            next_bigH += fraction * hc
        return next_bigK, next_bigH, (next_bigK / next_bigH)

    def find_nextk(self, krange, kstep):
        """
        Calculates aggregates according to PFE rules.
        calculates a projected littleK for each input littlek.
        Return k that gives smallest absolute difference.
        """
        assert isinstance(kstep, float)
        diff = float("inf")
        best_k = float("-inf")
        for dummy_k in range(1,  krange * int(kstep)):
            k = 1.0 * dummy_k / kstep
            new_littleK = self.calc_littlek(k)[2]
            new_diff = abs(k - new_littleK)
            if new_diff < diff:
                diff = new_diff
                best_k = k

        assert diff < .01
        assert best_k != 1 / kstep, "min k"
        assert best_k != (1.0 * krange * int(kstep) - 1) / kstep, "max k"
        return best_k, diff

    def start_next(self, krange, kstep):
        """
        initializes next economy.
        """
        #PFE Solutions
        pfe_soln = self.find_nextk(krange, kstep)[0]
        new_K, new_H, new_k = self.calc_littlek(pfe_soln)
        assert abs(new_k - pfe_soln) < .01

        # Initialize Individuals
        poor_bequest, poor_physical, poor_hc, poor_edu = self.bequest_division(self._indiv_p, new_k)
        new_poor = Individual(self._indiv_p.get_beta(), self._indiv_p.get_theta_bar(),
                              poor_bequest, poor_edu, "p")
        rich_bequest, rich_physical, rich_hc, rich_edu =  self.bequest_division(self._indiv_r, new_k)
        new_rich = Individual(self._indiv_r.get_beta(), self._indiv_r.get_theta_bar(),
                              rich_bequest, rich_edu, "r")

        #Initialize Firm
        new_firm = CESProd(sigma=self._firm.get_sigma(), epsilon=self._firm.get_epsilon(),
                           K=new_K, H=new_H, A=self._firm._A)

        #Initialize_school
        new_school = self._school.copy()

        #Initialize Next Economy
        new_econ = Economy(new_poor, new_rich, new_school, new_firm, self._fraction_r, self._imperfections)
        return new_econ



    def find_khat(self):
        """
        Find the level of k such that an individual leaves a bequest even when b=0.
        """
        return self._firm.inverse_wages(self._indiv_p.calc_theta())

    def regimeI_find_bR(self, in_k):
        """
        A helper function that finds bR for any given k in regime I.
        This is necessary to find the SS in regime 1.
        An input to regimeI_findnextk.
        """
        out_e = self.optimal_e(in_k)
        tot_H = 1.0 * self._fraction_r * self._school.hc(out_e) + (1 - self._fraction_r)
        first = (tot_H * in_k) / self._fraction_r
        out_BR = first - out_e
        tot_K = self._fraction_r * (out_BR - out_e)
        return tot_K, tot_H, out_BR, out_e

    def regimeI_findnextk(self, k, krange, kstep):
        """
        A helper function that performs find_nextk assuming that the economy stays in regime I.
        An input to find_ss_regimeI.
        """
        new_K, new_H, new_bR, new_e = self.regimeI_find_bR(k)
        new_poor = Individual(bequest=0, theta_bar=self._indiv_p.get_theta_bar(), educ=0,
                                  beta=self._indiv_p.get_beta(), type = "p")
        new_rich = Individual(bequest=new_bR, theta_bar=self._indiv_r.get_theta_bar(), educ=new_e,
                                  beta=self._indiv_r.get_beta(), type = "r")
        new_firm = CESProd(sigma = self._firm.get_sigma(), epsilon=self._firm.get_epsilon(),
                               K=new_K, H=new_H, A=self._firm.get_A())
        new_school = self._school.copy()
        new_economy = Economy(indiv_p=new_poor, indiv_r=new_rich, school=new_school,
                                  firm=new_firm, fraction_r=self._fraction_r)

        #Next calc
        out_k = new_economy.find_nextk(krange=krange, kstep=kstep)[0]
        return out_k

    def find_ss_regimeI(self, krange, kstep):
        """
        Calculates the states state in regime 1.
        Steps: guess k --> find bR -- > find nextk.
        Return the minimum, which is the steady state.
        """
        assert isinstance(kstep, float)
        diff = float("inf")
        best_k = float("-inf")
        for dummy_k in range(1,  krange * int(kstep)):
            in_k = 1.0 * dummy_k / kstep
            newk = self.regimeI_findnextk(k=in_k, krange=krange, kstep= kstep)
            new_diff = abs(newk - in_k)
            if new_diff < diff:
                diff = new_diff
                best_k = in_k
        assert diff < .01
        assert best_k != 1 / kstep, "min k"
        assert best_k != krange / kstep, "max k"
        return best_k, diff

    def no_regimeI_ss(self, krange, kstep):
        """
        return true if khat > regime I ss level of k.
        """
        return self.find_ss_regimeI(krange, kstep) >= self.find_khat()

    def avoid_trivial_ss(self):
        """
        An assertion that the economy avoids the trivial steady state.
        """
        pass

    def regimeIII_ss(self, krange, kstep, ktol):
        """
        Find steady state in regime 3.
        """
        pass
