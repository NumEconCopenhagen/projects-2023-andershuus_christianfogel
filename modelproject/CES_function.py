from scipy import optimize
import numpy as np
from types import SimpleNamespace
from scipy.optimize import minimize



class NumericalSolutionCES():

    def __init__(self):
        """
        parameters of the model:
        sigma: The degree of substitutability between consumption and leisure (Cannot be 1)
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

        par.beta = 0.5
        par.A = 20
        par.L = 75
        par.sigma=0.99

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
        sol_h = optimize.minimize(self.firm_profit,x0,args = (p,),bounds=bound,method='L-BFGS-B')

        #c. unpack solution
        sol.h_star = sol_h.x[0]
        sol.y_star = self.production_function(sol.h_star)
        sol.pi_star = p*sol.y_star-sol.h_star

        #d. Save the results. 
        return sol.h_star, sol.y_star, sol.pi_star

    def utility(self,c,l):

        """utility of the consumer (CES utility function)
        
        Args:
        c
        l
        sigma

        Returns:
        u: utility of consumer()
        
        """

        #a. unpack
        par = self.par

        #b. output
        return -(c**((par.sigma-1)/par.sigma)+l**((par.sigma-1)/par.sigma))**(par.sigma/(par.sigma-1))



    def utility_optimize(self,x):

        """utility_optimize

        Args:
        p

        Calls functions:
        utility()

        Returns:
        utility function
        """
        return self.utility(x[0],x[1])

    def ineq_constraint(self,x,p):

        """ineq_constraint

        Args:
        p
        x

        Calls functions:
        profit_maximization()

        Returns:
        sol.h_star: Optimal working hours for a given p
        Constraint for the income. The income must be positive, which is the interpretation of the return. 
        """
        #a. unpack
        par = self.par

        #We are intersted in the profit (pi) for a given price. 
        pi_constraint = self.firm_profit_maximization(p)[2]

        #This must be higher or equal to 0. 
        return pi_constraint+par.L-(p*x[0]+x[1]) # violated if negative

    def maximize_utility(self,p): 

        """maximize_utility

        Args:
        p
        L


        ineq_constraint(): Constraint for the income. The income must be positive.
        utility_optimize(): utility function

        Returns:
        c_star: Optimal consumption for a given p
        l_star: Optimal leisure for a given p
        """

        par = self.par

        # a. setup

        #The bounds as before. Consumption cannot be negative. 
        bounds = ((0,np.inf),(0,par.L))
        ineq_con = {'type': 'ineq', 'fun': self.ineq_constraint,'args': (p,)} 

        # b. call optimizer
        x0 = (25,8) # fit the equality constraint

        #Important that we have the constraint as this is the condition in the problem. 
        result = minimize(self.utility_optimize,x0,

                                    method='SLSQP',

                                    bounds=bounds,

                                    constraints=[ineq_con],

                                    options={'disp':False})

        c_star, l_star = result.x

        #c. output
        return c_star, l_star

    def market_clearing(self,p):

        """market_clearing

        Args:
        p
        L

        firm_profit_maximization
        maximize_utility

        Returns:
        goods_market_clearing: Market clearing for goods
        labor_market_clearing: Market clearing for labor
        """

        #a. unpack
        par = self.par

        #b. optimal behavior of firm for a given price
        h,y,pi=self.firm_profit_maximization(p)

        #c. optimal behavior of consumer for a price
        c,l=self.maximize_utility(p)

        #d. market clearing
        goods_market_clearing = y - c
        labor_market_clearing = h - par.L + l


        #e. output
        return goods_market_clearing, labor_market_clearing

    def find_relative_price(self,p_lower=0,p_upper=100):

        """find_relative_price (in market clearing)

        Args: 
        p_lower: lower bound for relative price
        p_upper: upper bound for relative price

        function_to_solve
        self.market_clearing

        Returns:
        p: relative price in market clearing
        good_clearing: Should be 0 in order to have clearing
        labor_clearing: Should be 0 in order to have clearing
        """

        # a. unpack
        par = self.par

        # b. Define the clearing function, which we want to be 0. We use Walras' law as we then only need to clear one market. 
        def function_to_solve(p):
            return self.market_clearing(p)[0]

        # c. We are finding when the market is clearing. brentq has found be faster for this case compared to bisect. 
        result = optimize.root_scalar(function_to_solve, method='brentq', bracket=[p_lower, p_upper])

        # d. check if a solution was found. We do only want to save the results in this case. 
        if result.converged:
            #Save the result from the root finding. 
            p = result.root
            good_clearing=self.market_clearing(p)[0]
            labor_clearing=self.market_clearing(p)[1]
            print(f' Beta = {par.beta:.2f}. Sigma = {par.sigma:.2f}  p = {p:.2f} -> Good clearing = {good_clearing:.2f}. Labor clearing = {labor_clearing:.2f}.')
        else:
            print("Fail")
            p, good_clearing, labor_clearing = None, None, None

        return p, good_clearing, labor_clearing