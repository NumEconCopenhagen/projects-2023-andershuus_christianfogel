from scipy import optimize
import numpy as np
from types import SimpleNamespace
from scipy.optimize import minimize



class NumericalSolutionCES():

    def __init__(self):

        "Initial parameters. They can be changed as seen with beta and sigma"

        par = self.par = SimpleNamespace()

        # a. parameters

        par.beta = 0.5
        par.A = 20
        par.L = 75
        par.sigma=0.99

        # b. solution

        sol = self.sol = SimpleNamespace()

    def production_function(self,h):

        "production function of the firm. Specified in the text"

        #a. unpack

        par = self.par

        #b. production function

        y = par.A*h**par.beta

        #c. output

        return y

    def firm_profit(self,h,p):

        "profit function of the firm. The wage is normalized to 1"

        #a. profit
        pi = p*self.production_function(h)-h

        #b. output
        return -pi

    def firm_profit_maximization(self,p):

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

        "CES-utility function"

        #a. unpack
        par = self.par

        #b. output
        return -(c**((par.sigma-1)/par.sigma)+l**((par.sigma-1)/par.sigma))**(par.sigma/(par.sigma-1))



    def utility_optimize(self,x):

        #Rewriting the function in order to get a constraint later on. 
        par = self.par
        return self.utility(x[0],x[1])

    def ineq_constraint(self,x,p):

        #a. unpack
        par = self.par

        #We are intersting in the profit (pi) for a given price. 
        h_constraint, y_constraint, pi_constraint = self.firm_profit_maximization(p)

        #This must be higher or equal to 0. 
        return pi_constraint+par.L-(p*x[0]+x[1]) # violated if negative

    def maximize_utility(self,p): 

        "maximize utility using an optimizer"

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

        "calculating the excess demand of the good and working hours"

        #a. unpack
        par = self.par
        sol = self.sol

        #b. optimal behavior of firm for a given price
        h,y,pi=self.firm_profit_maximization(p)

        #c. optimal behavior of consumer for a price
        c,l=self.maximize_utility(p)

        #d. market clearing
        goods_market_clearing = y - c
        labor_market_clearing = h - par.L + l


        #e. output
        return goods_market_clearing, labor_market_clearing

    def find_relative_price(self,tol=1e-4,iterations=500, p_lower=0.25, p_upper=0.75):

        "find price that causes markets to clear. We use Walras law and therefore focus on the good market (Same result to do for labor market)"

        # a. unpack
        par = self.par
        sol = self.sol

        #Initial values.                                                                                                       
        i=0

        #Want a max for the iterations. The whole loop is based on the algorithm. 

        while i<iterations:

            #Define the mean price. 

            p=(p_lower+p_upper)/2

            #Find the function value. 

            f = self.market_clearing(p)[0]

            #Criteria for stopping. 
            if np.abs(f)<tol: 
                good_clearing=self.market_clearing(p)[0]
                labor_clearing=self.market_clearing(p)[1]
                consumption=self.maximize_utility(p)[0]
                print(f' Step {i:.2f}: Beta = {par.beta:.2f}. Sigma = {par.sigma:.2f}  p = {p:.2f} -> Good clearing = {good_clearing:.2f}. Labor clearing = {labor_clearing:.2f}. Consumption = {consumption:.2f}')
                break

            #First option in step 6 for the algortithm. 
            elif self.market_clearing(p_lower)[0]*f<0:
                p_upper=p

            #Second option in step 6. 
            elif self.market_clearing(p_upper)[0]*f<0:

                p_lower=p

            #To stop if the loop is broken. 

            else: 
                print("Fail")
                return None

            #Updating i, by adding 1. 
            i+=1

        return p, good_clearing, labor_clearing