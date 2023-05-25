# Importing relevant packages
import sympy as sm
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"-"}) 
from scipy import optimize
from types import SimpleNamespace


class Question_1_4:
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

    def solve_social_optimum(self):
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
        return L_optimal, tau_optimal, utility_optimal

class Question_1_5():
    def __init__(self,sigma,rho):
        self.par = SimpleNamespace()
        
        # Parameters
        self.par.alpha = 0.5
        self.par.kappa = 1.0
        self.par.nu = (2.0 * 16.0 ** 2.0)**(-1.0)
        self.par.w = 1.0
        self.par.epsilon=1
        self.par.tau=0.52261 #Value from question 1.4 with 5 decimals. 

        self.par.sigma=sigma
        self.par.rho=rho

    def consumption_constraint(self, L):
        par = self.par
        return par.kappa + (1 - par.tau) * par.w * L
    
    def disutility(self,L): 
         par = self.par
         return par.nu*L**(1+par.epsilon)/(1+par.epsilon)

    def utility(self, G,L):
        par = self.par

        power1 = (par.sigma-1)/par.sigma
        power2= par.sigma/(par.sigma-1)

        numerator=((par.alpha*self.consumption_constraint(L)**power1+(1-par.alpha)*G**power1)**power2)**(1-par.rho)-1

        return -numerator/(1-par.rho)+self.disutility(L)

    def maximize_utility(self, G):
            result = optimize.minimize_scalar(lambda L: self.utility(G, L), bounds=(0, 24), method='bounded')
            optimal_L = result.x
            max_utility = -result.fun
            return optimal_L
    
    def government_constraint(self, L):
        par = self.par
        return par.tau * par.w * L


    def clearing(self,G):

        #a. unpack
        par = self.par

        #b. optimal behavior of firm for a given price
        L_opt=self.maximize_utility(G)

        #c. optimal behavior of consumer for a price
        G_clearing=self.government_constraint(L_opt)

        #d. market clearing
        clearing = G-G_clearing

        #e. output
        return clearing
    
    def find_government_consumption(self,G_lower=1,G_upper=100):
        # a. unpack
            par = self.par

            # b. Define the clearing function, which we want to be 0. Step 2 in described algorithm. 
            def function_to_solve(G):
                return self.clearing(G)

            # c. Step 3 in described algorithm.
            result = optimize.root_scalar(function_to_solve, method='brentq', bracket=[G_lower, G_upper])

            # d. Step 4 in described algortithm. 
            if result.converged:
                #Save the result from the root finding. 
                G_solution = result.root
                clearing=self.clearing(G_solution)
                L_optimal=self.maximize_utility(G_solution)
                print(f' sigma = {par.sigma:.3f}. rho = {par.rho:.3f}  G = {G_solution:.6f} -> clearing = {clearing:.2f}.Working hours = {L_optimal:.2f}')
            else:
                print("Fail")
                G_solution, clearing = None, None

            #return G_solution, clearing


class Question_1_6():
    def __init__(self,sigma,rho):
        self.par = SimpleNamespace()

        # Parameters
        self.par.alpha = 0.5
        self.par.kappa = 1.0
        self.par.nu = (2.0 * 16.0 ** 2.0)**(-1.0)
        self.par.w = 1.0
        self.par.epsilon=1
        self.par.sigma=sigma
        self.par.rho=rho
    def consumption_constraint(self,L,tau):

        par = self.par
        
        return par.kappa+(1-tau)*par.w*L

    def government_constraint(self, L, tau):
        par = self.par
        return tau * par.w * L

    def disutility(self,L): 
         par = self.par
         return par.nu*L**(1+par.epsilon)/(1+par.epsilon)

    def target(self, params):
        par = self.par
        L_optimal,tau_optimal = params

        #Running the regression with the target values, where the target values are those to be optimized. 
        C_optimal = self.consumption_constraint(L_optimal,tau_optimal)
        G_optimal =  self.government_constraint(L_optimal,tau_optimal)
        disutility_optimal = self.disutility(L_optimal)

        power1 = (par.sigma-1)/par.sigma
        power2= par.sigma/(par.sigma-1)

        numerator=((par.alpha*C_optimal**power1+(1-par.alpha)*G_optimal**power1)**power2)**(1-par.rho)-1

        return -numerator/(1-par.rho)+disutility_optimal
    def solve_social_optimum_Q6(self):
        """ This is to find optimal value of alpha and sigma for question 6
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
        #return L_optimal, tau_optimal, utility_optimal