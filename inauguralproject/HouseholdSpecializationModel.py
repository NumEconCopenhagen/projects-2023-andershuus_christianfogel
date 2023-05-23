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

        # d. disutlity of work
        epsilon_ = 1+1/par.epsilon
        TM = LM+HM
        TF = LF+HF
        disutility = par.nu*(TM**epsilon_/epsilon_+TF**epsilon_/epsilon_)
        
        return utility - disutility

    def solve_discrete(self):
        """ This code solves the model distretly
        See utility for code input. We only assume that the working hours can be divided
        in half hours intervals. 
        
        """
        
        par = self.par
        opt = self.opt
        
        # a. all possible choices as defined above. 
        x = np.linspace(0+1e-8,24,49)
        #note, we start from a small number just above 0, because when sigma = 0.5 then the exponent of H
        # will be equal to -1, and we cannot have 1/0.

        LM,HM,LF,HF = np.meshgrid(x,x,x,x) # all combinations in one meshgrid. 
    
        # Making the vectors for the 4 variables. 
        LM = LM.ravel() 
        HM = HM.ravel()
        LF = LF.ravel()
        HF = HF.ravel()

        # b. calculate utility for all possible outcomes. 
        u = self.calc_utility(LM,HM,LF,HF)
    
        # c. set to minus infinity if constraint is broken. Cannot work more than 24 hours per day
        I = (LM+HM > 24) | (LF+HF > 24) # | is "or"
        u[I] = -np.inf
    
        # d. finding index for maximizing argument
        j = np.argmax(u)
        
        #Finding maximizing values
        opt.LM = LM[j]
        opt.HM = HM[j]
        opt.LF = LF[j]
        opt.HF = HF[j]

        #Calculate ratio between hours worked at home for men and women
        opt.HF_HM = HF[j]/HM[j] 

        return opt
    
    def print_table_q1(self):
        """ This code prints a cross tabel for the discrete solution given in question 2 for different values of alpha and sigma. 
        See utility for code input. 
        """

        # a. load parameters
        par = self.par

        # b. empty text for a start
        text = ''
        
        # c. write alpha/Sigma in the corner of the code
        text += f'{"Alpha/Sigma":<7s}{"":1s}'

        # b. making the sigma values as the column values. 
        for sigma in np.linspace(0.5,1.5,3):
            par.sigma=sigma
            text += f'{sigma:8.2f}'
        text += '\n' + '-'*40 + '\n' # we add horizontal separator
        
        # d. Making the rows. 
        for i,alpha in enumerate(np.linspace(0.25,0.75,3)):
            par.alpha=alpha 
            
            #We need i>0 to avoid a blank row as the first row, so a layout fix. 
            if i > 0: 
                text += '\n'
            
            #alpha as row-values and with a vertical separator
            text += f'{alpha:10.2f} |' 

            #Plotting values of HF/HM
            for sigma in np.linspace(0.5,1.5,3): 
                par.sigma=sigma
                sol = self.solve_discrete() #call the solve function
                text += f'{sol.HF_HM:8.2f}' #plot values 
        
        #e. reset values of alpha and sigma. This is necessary to avoid problems later in the following questions. 
        par.alpha = 0.5
        par.sigma = 1.0
        
        # f. Printing the table
        print(f"Table of HF/HM values:\n{text}")

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
        opt.HF_HM = solution_cont.x[3]/solution_cont.x[1] 

        return opt
    
    def solve_wF_vec(self, discrete=False):
        """Solve model for vector of female relative wages to men as men's wage is numeraire. 
        For discrete case we use discrete=True and discrete=False for the continuous model. 
        This is also step 2 in the described algorithm for question 4. 
        """

        par = self.par
        opt = self.opt

        # We create a vector to store the log ratio of HF/HM for the different log ratios of wF/wM
        log_HF_HM = np.zeros(par.wF_vec.size)

        # We loop over each value of wF in wF_vec
        for i, wF in enumerate(par.wF_vec):
            par.wF = wF # Set the new value of wF
            
            # Solve the model based on whether we have a discrete or continuous case
            if discrete==True:
                opt = self.solve_discrete()
                log_HF_HM[i] = np.log(opt.HF_HM)
            else:
                opt = self.solve_continuously()
                log_HF_HM[i] = np.log(opt.HF_HM)

        #We return wF to original value
        par.wF = 1.0
        return log_HF_HM #Return the vector of log ratio of HF/HM


    #Defining the target function for question 4. Step 1 in the algorithm described in section 4. 
    def target(self, params):
        par = self.par
        alpha_target, sigma_target = params
        beta0_target=0.4
        beta1_target=-0.1

        #Running the regression with the target values, where the target values are those to be optimized. 
        beta0,beta1=self.run_regression(alpha_target,sigma_target)
        return (beta0_target-beta0)**2+(beta1_target-beta1)**2
    
    #Step 3 in described algorithm. 
    def run_regression(self,alpha_optimal, sigma_optimal):
        """ This regression is for question 4 and 5 
        Our input is the value from the solve_wF_Vec
        Therefore it can be used for discrete and continuous case"""

        par = self.par
        opt = self.opt

        par.alpha=alpha_optimal
        par.sigma=sigma_optimal

        #Defining the x and y values for the regression.
        x = np.log(par.wF_vec)
        y = self.solve_wF_vec(discrete=False)
        A = np.vstack([np.ones(x.size),x]).T

        #Making the regression and returning the parameter estimates. 
        opt.beta0,opt.beta1 = np.linalg.lstsq(A,y,rcond=None)[0]

        return opt.beta0,opt.beta1

    #Step 4 in the described algoirthm for question 4. 
    def solve_optimal_alpha_sigma(self,used_seed):
        """ This is to find optimal value of alpha and sigma for question 4
        We use Nelder-Mead because we do not have constraints"""
        #Setting a random seed. We will loop over the seed in order to make robustness checks. 
        np.random.seed(used_seed)

        #Defining the initial guess based on the seed. 
        init_guess = [np.random.uniform(0.01,1), np.random.uniform(0.01,1)]

        #Bounds for alpha and sigma. 
        bounds = [(0.01, 1-1e-8), (0.01, 100)]

        result = optimize.minimize(self.target,init_guess, method='Nelder-Mead', bounds=bounds)

        #Saving result from the optimizer. 
        alpha_result=result.x[0]
        sigma_result=result.x[1]

        #Saving the results in one parameter in order to include in target function. 
        params_result=alpha_result,sigma_result

        #Returning the values. 
        return alpha_result, sigma_result, self.target(params_result)   

