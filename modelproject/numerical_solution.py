from scipy import optimize
import numpy as np
from types import SimpleNamespace

class NumericalSolution():
    
    def __init__(self):
        """
        parameters of the model:
        alpha: Relative preferences towards consumption
        beta: Returns to scale
        A: Total factor productivity
        L: Labor endowment

        variables of the model:
        h: working hours
        l: leisure
        c: consumption
        u: utility
        y: total production
        p: relative price
        pi: firm profits

        """
        # a. create SimpleNamespaces
        par = self.par = SimpleNamespace()
        sol = self.sol = SimpleNamespace()

        # b. parameters
        par.alpha = 0.5
        par.beta = 0.5
        par.A = 20
        par.L = 75

    def production_function(self,h):
        """production function of the firm
        
        Args:
        h

        Returns:
        y: total production

        """

        #a. unpack
        par = self.par

        #b. production function
        y = par.A*h**par.beta
        
        #c. output
        return y

    def firm_profit(self,h,p):
        """profit function of the firm
        
        Args:
        h
        p
        
        Calls function:
        production_function()

        Returns:
        -pi: negative value of firm profits

        The reason for return a negative value of firm profits is so we can use a minimizer
        for firm profit maximization.
        
        """

        #a. profit
        pi = p*self.production_function(h)-h

        #b. output
        return -pi
    
    def firm_profit_maximization(self,p):

        """firm profit maximization

        Args:
        p

        Calls functions:
        firm_profit()
        production_function()

        Returns:
        sol.h_star: Optimal working hours for a given p
        sol.y_star: Optimal total production for a given p
        sol.pi_star: Firm profits for optimal production and working hours and for a given p
        
        """

        #a. unpack
        par = self.par
        sol = self.sol

        #b. call optimizer. h cannot be higher than the initial endowment
        bound = ((0,par.L),)
        x0=[0.0]
        #Other optimizer could also have been used. 
        sol_h = optimize.minimize(self.firm_profit,x0,args = (p,),bounds=bound,method='L-BFGS-B')

        #c. unpack solution
        sol.h_star = sol_h.x[0]
        sol.y_star = self.production_function(sol.h_star)
        sol.pi_star = p*sol.y_star-sol.h_star

        return sol.h_star, sol.y_star, sol.pi_star

    def utility(self,c,h):
        """utility of the consumer
        
        Args:
        c
        h
        alpha

        Returns:
        u: utility of consumer
        
        """

        #a. unpack
        par = self.par

        #b. utility (cobb-douglas)
        u = c**par.alpha*(par.L-h)**(1-par.alpha)

        #c. output
        return u

    def income(self,p):
        """consumer's income/budget constraint
        
        Args:
        p

        Calls function:
        firm_profit_maximization()

        Returns:
        sol.Inc: Consumers income for a given p

        """

        #a. unpack
        par = self.par
        sol = self.sol

        #b. budget constraint. Minus because the income is defined negatively in order to optimize. 
        pi_inc=self.firm_profit_maximization(p)[2]
        sol.Inc = pi_inc+par.L

        #c. output
        return sol.Inc

    def maximize_utility(self,p):

        """utility maximization of the consumer
        
        Args: p

        Calls function:
        income()

        Returns:
        sol.c_star: optimal consumption for a given p
        sol.l_star: optimal leisure for a given p
        
        """
        
        # a.unpack
        par = self.par
        sol = self.sol

        # a. solve using standard solutions
        utility_inc=self.income(p)

        sol.c_star = par.alpha*utility_inc/p
        sol.l_star = (1-par.alpha)*utility_inc

        return sol.c_star, sol.l_star
    
    def market_clearing(self,p):
        """calculating excess supply of the good and excess demand of labor

        Args:
        p

        Functions called:
        firm_profit_maximization()
        maximize_utility()

        Returns:
        goods_market_clearing: excess supply of the good for a given p
        labor_market_clearing: excess demand of labor for a given p
        
        """
        #a. unpack
        par = self.par

        #b. optimal behavior of firm
        h,y,pi=self.firm_profit_maximization(p)

        #c. optimal behavior of consumer
        c,l=self.maximize_utility(p)

        #b. market clearing
        goods_market_clearing = y - c
        labor_market_clearing = h - par.L + l

        return goods_market_clearing, labor_market_clearing
    
    def find_relative_price(self,tol=1e-4,iterations=500, p_lower=1e-4, p_upper=100):
        """finding the optimal relative price

        Args: Step 2 in the algorithm is specifying this. 
        tol: tolerance level for market clearing. Default value set to 1e-4
        iterations: number of iterations allowed for while loop. Default value set to 500
        p_lower: lowest bound for optimal relative prices. Default value set to 1e-4. 
        -We set it to a small number above 0 to avoid deviding by 0 in the maximize_utility() method.
        p_upper: upper bound for optimal relative prices. Default set to 100.
        -We set it to 100, just to have a large number as our upper bound.

        Functions called:
        market_clearing()

        Returns:
        p: optimal relative price
        good_clearing: value of excess supply at the optimal relative price
        labor_clearing: value of excess demand for labor at the optimal relative price

        """

        #Initial values.                                                                                                       
        i=0

        #step 3-7 in algorithm
        while i<iterations: 
            
            #step 4 in algorithm
            p=(p_lower+p_upper)/2

            #step 5 in algorithm
            f = self.market_clearing(p)[0]

            #step 3+6-7 in algorithm
            if np.abs(f)<tol: 
                good_clearing=self.market_clearing(p)[0]
                labor_clearing=self.market_clearing(p)[1]
                print(f' Step {i:.2f}: p = {p:.2f} -> Good clearing = {good_clearing:.8f}. Labor clearing = {labor_clearing:.2f}. ')
                break
            elif self.market_clearing(p_lower)[0]*f<0:
                p_upper=p
                #print(f' Step {i:.2f}: p = {p:.2f} -> {f:12.8f}')
            elif self.market_clearing(p_upper)[0]*f<0:
                p_lower=p
                #print(f' Step {i:.2f}: p = {p:.2f} -> {f:12.8f}')
            else: 
                print("Fail")
                return None
            
            i+=1
        return p, good_clearing, labor_clearing