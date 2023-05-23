from types import SimpleNamespace

import numpy as np
from scipy import optimize


class HouseholdSpecializationModelClass:
    
    def __init__(self):
        """ We are having the following variables
        rho
        nu 
        epsilon
        omega
        psi - measure relative disutility for men working at home. 
        alpha - Decides the relative productivity of men and women in the home
        sigma - Decides the degree of substitution. 
        wM - Wages for men. A numeraire
        wF - Wages for women 
        wF_vec - Vector for women's wage used for question 2 and onwaard
        beta0_target - From the paper for question 4
        beta1_target - From the paper for question 4
          
            """

        # a. create namespaces
        par = self.par = SimpleNamespace()
        opt = self.opt = SimpleNamespace()

        # b. preferences
        par.rho = 2.0
        par.nu = 0.001
        par.epsilon = 1.0
        par.omega = 0.5 
        par.psi = 1.0

        # c. household production
        par.alpha = 0.5
        par.sigma = 1.0

        # d. wages
        par.wM = 1.0
        par.wF = 1.0
        par.wF_vec = np.linspace(0.8,1.2,5)
    
        # e. targets
        par.beta0_target = 0.4
        par.beta1_target = -0.1

    def calc_utility(self,LM,HM,LF,HF):
        """ This function is calculating the utility for the household
        The input variables are: 
        wM
        LM - Hours at work for men
        HM - Hours at home for men
        wF
        LF - Hours at work for women
        HF - Hours at home for women
        Sigma
        Alpha
        H = Consumption of home production
        Omega
        Rho
        Q = Total consumption
        Epsilon
        TM - Total working hours for men
        TF - Total working hours for women
        Psi
        Nu

        """

        par = self.par

        # a. consumption of market goods
        C = par.wM*LM + par.wF*LF

        # b. home production 
        #Based on the value of sigma we use the specified function for H
        if par.sigma == 0:
            H=np.fmin(HM,HF)
        elif par.sigma == 1:
            H = HM**(1-par.alpha)*HF**par.alpha
        else:
            H= ((1-par.alpha)*HM**((par.sigma-1)/par.sigma)+par.alpha*HF**((par.sigma-1)/par.sigma))**(par.sigma/(par.sigma-1))

        # c. total consumption utility
        Q = C**par.omega*H**(1-par.omega)
        utility = np.fmax(Q,1e-8)**(1-par.rho)/(1-par.rho)

        # d. disutility of work. This is where we add our new parameter psi.
        epsilon_ = 1+1/par.epsilon
        TM = LM+HM
        TF = LF+HF
        disutility = par.nu*(TM**epsilon_/epsilon_+TF**epsilon_/epsilon_)+par.psi*HM
        
        return utility - disutility

    #We need negative value of utility function for the continuously function. 
    """Calculate the negative value of the utility function. Inputs are the same as in calc_utility()
    We multiply the utility by 100 to avoid imprecise estimates.
    """
    def utility_function(self, L): 
        return -self.calc_utility(L[0],L[1],L[2],L[3])*100
    
    def solve_continuously(self):
        """
        This code solves the model with continuously choice variables for working hours (at work and home)
        We use SLSQP for optimizing because we want to insert bounds and constraints.
        Input are the same as for the discrete case"""
        par = self.par
        opt = self.opt
        
        #Define the bounds and constraints. 
        constraint_men = ({'type': 'ineq', 'fun': lambda L:  24-L[0]-L[1]})
        constraint_women = ({'type': 'ineq', 'fun': lambda L:  24-L[2]-L[3]})
        
        #Hours worked cannot be 0 as we will divide with 0 when sigma is equal to 0.5. Therefore, leisure can also just be below 24 hours.
        bounds=((0,24-1e-8),(0+1e-8,24), (0,24-1e-8), (0+1e-8,24))
        initial_guess=[6,6,6,6]

        # Call optimizer based on the input above. 
        solution_cont = optimize.minimize(
        self.utility_function, initial_guess,
        method='SLSQP', bounds=bounds, constraints=(constraint_men, constraint_women))
        
        # Save results. x[0]=L[0] and so on. 
        opt.LM = solution_cont.x[0]
        opt.HM = solution_cont.x[1]
        opt.LF = solution_cont.x[2]
        opt.HF = solution_cont.x[3]

        #Calculating the ratio
        opt.HF_HM = solution_cont.x[3]/solution_cont.x[1] #calculate ratio

        return opt
    
    def solve_wF_vec(self):
        """ Solve model for vector of female wages. The same as in HouseholdSpecializationModel except we do not have the
        discrete case here. 
        This is also step 2 in the described algorithm for question 4.
        """

        par = self.par
        opt = self.opt

        # We create a vector to store the log ratio of HF/HM for the different log ratios of wF/wM
        log_HF_HM = np.zeros(par.wF_vec.size)

        # We loop over each value of wF in wF_vec
        for i, wF in enumerate(par.wF_vec):
            par.wF = wF # Set the new value of wF
            
            # Solve the model
            opt = self.solve_continuously()
            log_HF_HM[i] = np.log(opt.HF_HM)

        par.wF = 1.0
        #We return wF to original value
        return log_HF_HM #Return the vector of log ratio of HF/HM


    #Defining the target function for question 5. Step 1 in the algorithm described in section 4. 
    def target(self, params):
        par = self.par
        sigma_target,psi_target = params
        beta0_target=0.4
        beta1_target=-0.1

        #Running the regression with the target values, where the target values are those to be optimized. 
        beta0,beta1=self.run_regression(sigma_target,psi_target)
        return (beta0_target-beta0)**2+(beta1_target-beta1)**2
    
    #Step 3 in described algorithm. 
    def run_regression(self, sigma_optimal,psi_optimal):
        """ This regression is for question 5 
        Our input is the value from the solve_wF_Vec"""

        par = self.par
        opt = self.opt

        par.sigma=sigma_optimal
        par.psi=psi_optimal


        x = np.log(par.wF_vec)
        y = self.solve_wF_vec()
        A = np.vstack([np.ones(x.size),x]).T

        #Making the regression and returning the parameter estimates. 
        opt.beta0,opt.beta1 = np.linalg.lstsq(A,y,rcond=None)[0]

        return opt.beta0,opt.beta1

    #Step 4 in the described algoirthm. 
    def solve_optimal_sigma_psi(self,used_seed):
        """ This is to find optimal value of sigma and psi for question 5
        We use Nelder-Mead because we do not have constraints"""
        #Setting a random seed. We will loop over the seed in order to make robustness checks. 
        np.random.seed(used_seed)

        #Defining the initial guess based on the seed. 
        init_guess = [np.random.uniform(0.01,1)]

        #Bounds for sigma and psi. 
        bounds = [(0.01, 100), (-10,10)]

        #Optimize. 
        result = optimize.minimize(self.target,init_guess, method='Nelder-Mead', bounds=bounds)

        #Saving result from the optimizer. 
        sigma_result=result.x[0]
        psi_result=result.x[1]

        #Saving the results in one parameter in order to include in target function. 
        params_result=sigma_result,psi_result

        #Returning the values. 
        return sigma_result,psi_result, self.target(params_result)   

