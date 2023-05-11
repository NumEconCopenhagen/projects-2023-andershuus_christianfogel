import platformdirs
from scipy import optimize
import numpy as np
import ipywidgets as widgets # Interactive plots
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy.optimize import minimize

class numerical_solution_ces():

    def __init__(self):

        par = self.par = SimpleNamespace()

        # a. parameters
        par.alpha = 0.5
        par.beta = 0.5
        par.A = 20
        par.L = 75
        par.sigma=0.99

        # b. solution
        sol = self.sol = SimpleNamespace()
        
        sol.p = 1 # output price
        sol.w = 1 # wage
        sol.c = 1 # consumption
        sol.l = 1 # leisure
        sol.h = 1 # working hours
        sol.y = 1 # production
        sol.pi = 1 # profit
        sol.Inc = 1 # income

    def production_function(self,h):
        "production function of the firm"

        #a. unpack
        par = self.par

        #b. production function
        y = par.A*h**par.beta
        
        #c. output
        return y

    def firm_profit(self,h,p):
        "profit function of the firm"

        #a. profit
        pi = p*self.production_function(h)-h

        #b. output
        return -pi
    
    def firm_profit_maximization(self,p):

        #a. unpack
        par = self.par
        sol = self.sol

        #b. call optimizer
        bound = ((0,par.L),)
        x0=[0.0]
        sol_h = optimize.minimize(self.firm_profit,x0,args = (p,),bounds=bound,method='L-BFGS-B')

        #c. unpack solution
        sol.h_star = sol_h.x[0]
        sol.y_star = self.production_function(sol.h_star)
        sol.pi_star = p*sol.y_star-sol.h_star

        return sol.h_star, sol.y_star, sol.pi_star

    def utility(self,c,l):
        par = self.par
        return -(c**((par.sigma-1)/par.sigma)*l**((par.sigma-1)/par.sigma))**(par.sigma/(par.sigma-1))

    def utility_optimize(self,x):
        par = self.par
        return self.utility(x[0],x[1])
        
    def ineq_constraint(self,x,p):
        par = self.par
        h_constraint, y_constraint, pi_constraint = self.firm_profit_maximization(p)
        return pi_constraint+par.L-(p*x[0]+x[1]) # violated if negative
    def maximize_utility(self,p): 
        par = self.par
        # a. setup
        bounds = ((0,np.inf),(0,par.L))
        #constraint = ineq_constraint(p)
        ineq_con = {'type': 'ineq', 'fun': self.ineq_constraint,'args': (p,)} 

        # b. call optimizer
        x0 = (25,8) # fit the equality constraint
        result = minimize(self.utility_optimize,x0,
                                    method='SLSQP',
                                    bounds=bounds,
                                    constraints=[ineq_con],
                                    options={'disp':False})

        c_star, l_star = result.x
        
        return c_star, l_star
    

    def market_clearing(self,p):
        "calculating the excess demand of the good and working hours"
        #a. unpack
        par = self.par
        sol = self.sol

        #b. optimal behavior of firm
        h,y,pi=self.firm_profit_maximization(p)

        #c. optimal behavior of consumer
        c,l=self.maximize_utility(p)

        #b. market clearing
        goods_market_clearing = y - c
        labor_market_clearing = h - par.L + l

        return goods_market_clearing, labor_market_clearing
    
    def find_relative_price(self,tol=1e-4,iterations=500, p_lower=0.25, p_upper=0.75, delta=1e-8, adj=0.5):
        "find price that causes markets to clear"

        # a. unpack
        par = self.par
        sol = self.sol

        #Initial values.                                                                                                       
        i=0

        while i<iterations: 
            
            p=(p_lower+p_upper)/2
            f = self.market_clearing(p)[0]
            #Approximation of derivative. 
            #fp = (self.market_clearing(p+delta)[0]-f)/delta

            if np.abs(f)<tol: 
                good_clearing=self.market_clearing(p)[0]
                labor_clearing=self.market_clearing(p)[1]
                consumption=self.maximize_utility(p)[0]

                print(f' Step {i:.2f}: p = {p:.2f} -> {f:12.8f}. Good clearing = {good_clearing:.2f}. Labor clearing = {labor_clearing:.2f}. Consumption = {consumption:.2f}')
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