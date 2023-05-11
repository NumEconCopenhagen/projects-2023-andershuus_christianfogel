from scipy import optimize
import numpy as np
import ipywidgets as widgets # Interactive plots
import matplotlib.pyplot as plt
from types import SimpleNamespace


#Analytical solution, i.e. section 2
def consumption(alpha,beta,A, L=24): 
    """
        
    Args:

        alpha: relative preferences for consumption relative to leisure
        beta: output elasticity
        A: TFP
        
    Returns:
    
        p = relative price of consumption good (Wage is normalized to 1)
        h = working hours (labor demand)
        y = production
        pi = profit
        Inc = income for consumer
        c = consumption
        l = leisure 
        
    """
    # Initally find the price
    numerator=(L*alpha)**(1-beta)
    denominator = A*(beta**(beta/(1-beta))*(1-alpha)+alpha*beta**(1/(1-beta)))**(1-beta)
    p = numerator/denominator

    # Labor demand, production and profit for the firm
    h = (beta*p*A)**(1/(1-beta))
    y = A*(h)**beta
    pi = p*y - h

    # Define income, leisure and consumption
    Inc = pi+L
    c = alpha*Inc/p
    l = (1-alpha)*Inc

    #Check that good and labor market clear
    assert np.isclose(c, y,1e-8), 'Good market does not clear'
    assert np.isclose(L-h, l,0.0), 'Labor market does not clear'

    return p,h,y,pi,Inc,c,l

    #Print statement which are not used. 
    #print(f'price = {p:.2f}')
    #print(f'profit = {pi:.1f}')
    #print(f'income = {Inc:.1f}')
    #print(f'working = {h:.1f}')
    #print(f'leisure = {l:.1f}')
    #print(f'production = {y:.1f}')
    #print(f'consumption = {c:.1f}')

#Making an interactive plot of the solution
def interactive_figure(beta,A):
    """
        
    Args:

        beta: output elasticity
        A: TFP
        
    Returns:
    
        p_vec (in form of a figure)= vector of relative price of consumption good (Wage is normalized to 1) w.r.t. alpha
        c_vec (in form of a figure)= vector of consumption w.r.t. alpha
        l_vec (in form of a figure)= vector of leisure w.r.t. alpha
        
    """
    
    alpha_vec = np.linspace(1e-8,1-1e-8,10) #Cannot be 0 or 1 
    
    p_vec = np.empty(len(alpha_vec))
    h_vec = np.empty(len(alpha_vec))
    y_vec = np.empty(len(alpha_vec))
    pi_vec = np.empty(len(alpha_vec))
    Inc_vec = np.empty(len(alpha_vec))
    c_vec = np.empty(len(alpha_vec))
    l_vec = np.empty(len(alpha_vec))

    # a. Solving the model for the given alpha
    for i, alpha in enumerate(alpha_vec):
        p_vec[i],h_vec[i],y_vec[i],pi_vec[i],Inc_vec[i], c_vec[i], l_vec[i] = consumption(alpha,beta,A)
    
    # b. figure for price

    fig = plt.figure(figsize=(8, 10))
    ax = fig.add_subplot(3,1,1)
    ax.plot(alpha_vec, p_vec, label='Price of consumption relative to leisure')
    ax.set_xlim([0.05,0.95]) # 
    ax.set_ylim([0,20]) #
    ax.set_title("Price")
    ax.set_xlabel("Alpha") 
    ax.legend(loc= 'upper right')

    # c. figure for consumption
    bx = fig.add_subplot(3,1,2)
    bx.plot(alpha_vec, c_vec, label='Consumption')
    bx.set_xlim([0.05,0.95]) # 
    bx.set_ylim([0,70]) #
    bx.set_title("Consumption")
    bx.set_xlabel("Alpha") 
    bx.legend(loc= 'upper right')

    # d. figure for leisure
    cx = fig.add_subplot(3,1,3)
    cx.plot(alpha_vec, l_vec, label='Leisure')
    cx.set_xlim([0.05,0.95]) # 
    cx.set_ylim([0,70]) #
    cx.set_title("Leisure")
    cx.set_xlabel("Alpha") 
    cx.legend(loc= 'upper right')

    plt.tight_layout()

class numerical_solution():

    def __init__(self):

        par = self.par = SimpleNamespace()

        # a. parameters
        par.alpha = 0.5
        par.beta = 0.5
        par.A = 20
        par.L = 75

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

    def utility(self,c,h):
        "utility of the consumer"

        #a. unpack
        par = self.par

        #b. utility
        u = c**par.alpha*(par.L-h)**(1-par.alpha)

        #c. output
        return u

    def income(self):
        "consumer's income/budget constraint"

        #a. unpack
        par = self.par
        sol = self.sol

        #b. budget constraint. Minus because the income is defined negatively in order to optimize. 
        sol.Inc = sol.pi_star+par.L

        #c. output
        return sol.Inc


    def maximize_utility(self,p):
        
        # a.unpack
        par = self.par
        sol = self.sol

        # a. solve using standard solutions
        sol.c_star = par.alpha*sol.Inc/p
        sol.l_star = (1-par.alpha)*sol.Inc

        return sol.c_star, sol.l_star


    def utility_maximization(self,p): 

        #a. unpack
        par = self.par
        sol = self.sol

        #b. call optimizer
        #Bounds
        bounds = ((0,np.inf),(0,par.L))
        #Initial guess
        x0=[25,8]

        #Constraints. The income must be equal to or greater than the income. We first define l 
        constraint = sol.Inc-p*self.utility[0]-par.L-sol.h_star
        ineq_con = {'type': 'ineq', 'fun': constraint} 

    
        #sol_h = optimize.minimize(self.firm_profit,x0,args =(p,),bounds=bounds,method='L-BFGS-B')

        # b. call optimizer
        sol_con = optimize.minimize(self.utility,x0,
                             method='SLSQP',
                             bounds=bounds,
                             constraints=[ineq_con],
                             options={'disp':True})
        c_star = sol_con.x[0]
        l_star = sol_con.x[1]

        return c_star, l_star
    
    def market_clearing(self,p):
        "calculating the excess demand of the good and working hours"
        #a. unpack
        par = self.par
        sol = self.sol

        #b. optimal behavior of firm
        self.firm_profit_maximization(p)

        #c. optimal behavior of consumer
        self.maximize_utility(p)

        #b. market clearing
        goods_market_claering = sol.y_star - sol.c_star
        labor_market_clearing = sol.h_star - par.L + sol.l_star

        return goods_market_claering, labor_market_clearing
    
    def solve(self):
        "find price that causes markets to clear"

        # a. unpack
        par = self.par
        sol = self.sol

        #b. initial guess
        p_guess = 1.0

        #c. tolerance
        tolerance = 1e-8

        #d. iterations
        max_iterations=500
        i = 0

        





    
    def solve(self):

        #a. unpack
        par = self.par

        #b. initial guess
        p_guess = 1.0

        #c. tolerance
        tolerance = 1e-8

        #d. iterations
        max_iterations=500
        i = 0

        #e. solve for p


        





    def firm(self):
        "maximize firm profits"

        #a. unpack
        par = self.par
        sol = self.sol

        p = sol.p
        w = sol.w

        #b. objective function
        y = lambda h: par.alpha*h**par.beta

        obj = lambda h: -(p*y(h)-w*h)

        #c. call optimizer
        x0 = [0.0]
        res = optimize.minimize(obj,x0,bounds=((0,None),),method='L-BFGS-B')
        
        #d. unpack solution
        sol.h_star = res.x[0]
        sol.y_star = y(sol.h_star)


class lecture_5():
    
    def __init__(self):

        par = self.par = SimpleNamespace()

        # a. parameters
        par.alpha = 0.5
        par.beta = 0.5
        par.A = 20
        par.L = 24

        # b. grids
        par.grid_p = np.linspace(0.1,1.5,10)
        par.grid_mkt_clearing = np.zeros(10)

        # c. solution
        sol = self.sol = SimpleNamespace()
        
        sol.p = 1 # output price

    def firm(self):
        """ maximize firm profits """
            
        par = self.par
        sol = self.sol

        p = sol.p

        # a. solve
        f = lambda h: par.A*h**par.beta
        obj = lambda h: -(p*f(h)-h)
        x0 = [0.0]
        res = optimize.minimize(obj,x0,bounds=((0,par.L),),method='L-BFGS-B')
            
        # b. save
        sol.h_star = res.x[0]
        sol.y_star = f(sol.h_star)
        sol.pi = p*sol.y_star-sol.h_star

    def utility_c(self,c,h):
        """ utility of consumer """

        par = self.par

        return c**par.alpha*h**(1-par.alpha)

    def consumer(self):
        """ maximize utility of consumer """
        par = self.par
        sol = self.sol

        p = sol.p
        pi = sol.pi

        # a. solve
        obj = lambda h: -self.utility_c(par.alpha*(h+pi)/p,h) 
        res = optimize.minimize_scalar(obj,bounds=(0,1),method='bounded')
        
        # b. save
        sol.h_star = res.x
        sol.c_star = (sol.h_star+pi)/p
        sol.l_star = (1-par.alpha)*(pi+sol.h_star)

    def evaluate_equilibrium(self):
        """ evaluate equilirium """
        
        par = self.par
        sol = self.sol

        # a. optimal behavior of firm
        self.firm()

        # b. optimal behavior of households
        self.consumer()

        # c. market clearing
        sol.goods_mkt_clearing = sol.y_star - sol.c_star
        sol.labor_mkt_clearing = par.L - sol.l_star-sol.h_star
    
    def find_equilibrium(self):

        par = self.par
        sol = self.sol

        # a. grid search
        print('grid search:')
        for i,p in enumerate(par.grid_p):
            sol.p = p
            self.evaluate_equilibrium()
            par.grid_mkt_clearing[i] = sol.goods_mkt_clearing
            print(f' p = {p:.2f} -> {par.grid_mkt_clearing[i]:12.8f}')
        
        print('')

        # b. find bounds
        left = np.max(par.grid_p[par.grid_mkt_clearing < 0])
        right = np.min(par.grid_p[par.grid_mkt_clearing > 0])
        print(f'equilibrium price must be in [{left:.2f},{right:.2f}]\n')            

        # c. bisection search
        def obj(p):
            sol.p = p
            self.evaluate_equilibrium()
            return sol.goods_mkt_clearing

        res = optimize.root_scalar(obj,bracket=[left,right],method='bisect')
        sol.p = res.root
        print(f'the equilibrium wage is {sol.w:.4f}\n')

        # d. show result
        u = self.utility_c(sol.c_star,sol.h_star)
        print(f'capitalists  : c = {sol.c_star:6.4f}, h = {sol.h_star:6.4f}, u = {u:7.4f}')        
        print(f'goods market : {sol.goods_mkt_clearing:.8f}')
        print(f'labor market : {sol.labor_mkt_clearing:.8f}')