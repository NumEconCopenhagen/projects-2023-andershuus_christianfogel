# Importing relevant packages
import sympy as sm
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"-"}) 
from scipy import optimize
from types import SimpleNamespace


class WorkerClass:
    def __init__(self):
        self.par = SimpleNamespace()

        # Parameters
        self.par.alpha = 0.5
        self.par.kappa = 1.0
        self.par.nu = (2.0 * 16.0 ** 2.0)**(-1.0)
        self.par.w = 1.0
    
    def consumption_constraint(self,L,tau):

        par = self.par
        
        return par.kappa+(1-tau)*par.w*L

    def government_constraint(self, L, tau):
        par = self.par
        return tau * par.w * L

    def disutility(self,L): 
         par = self.par
         return par.nu*L**2/2

    def target(self, params):
        par = self.par
        L_optimal,tau_optimal = params

        #Running the regression with the target values, where the target values are those to be optimized. 
        C_optimal = self.consumption_constraint(L_optimal,tau_optimal)
        G_optimal =  self.government_constraint(L_optimal,tau_optimal)
        disutility_optimal = self.disutility(L_optimal)


        return -np.log(C_optimal**par.alpha*G_optimal**(1-par.alpha))+disutility_optimal

    def solve_optimal_alpha_sigma(self):
        """ This is to find optimal value of alpha and sigma for question 4
        We use Nelder-Mead because we do not have constraints"""

        #Defining the initial guess based on the seed. 
        init_guess = [(5,0.5)]

        #Bounds for alpha and sigma. 
        bounds = [(0.01, 24), (0.01,0.99)]

        result = optimize.minimize(self.target,init_guess, method='SLSQP', bounds=bounds)

        #Saving result from the optimizer. 
        L_optimal, tau_optimal = result.x

        #Saving the results in one parameter in order to include in target function. 
        params_result=L_optimal,tau_optimal

        #Finding the utility
        utility_optimal = -self.target(params_result)


        #Returning the values. 
        print(f' Optimal labor supply = {L_optimal:.2f}. Optimal tax rate = {tau_optimal:.2f}. Utility in optimum = {utility_optimal:.2f}.')
         