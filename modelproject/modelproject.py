from scipy import optimize
import numpy as np
import ipywidgets as widgets # Interactive plots
import matplotlib.pyplot as plt


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
        par.L = 24

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

    def firm(self,p,w):
        

        #a. unpack
        par = self.par

        #b. objective function
        y = lambda h: par.alpha*h**par.beta

        obj = lambda h: -(p*y-w*h)

        #c. call optimizer
        x0 = [0.0]
        res = optimize.minimize(obj,x0,bounds=((0,None),),method='L-BFGS-B') #look into method
    
        #d. unpack solution

        return res