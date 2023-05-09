from scipy import optimize
import numpy as np
import ipywidgets as widgets # Interactive plots
import matplotlib.pyplot as plt


#Analytical solution, i.e. section 2
def consumption(alpha,beta,A): 
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
    p = (24**(1-beta))/(beta*alpha)+((24*alpha)/(1-alpha))**(1-beta)*1/(A*beta**beta)

    # Labor demand, production and profit for the firm
    h = (beta*p*A)**(1/(1-beta))
    y = A*h**beta
    pi = p*y - h

    # Define income, leisure and consumption
    Inc = pi+24
    c = alpha*Inc/p
    l = (1-alpha)*Inc

    #Check that supply is equal to demand for both goods. Demand cannot be smaller due to assumption about monotonicity
    assert np.isclose(c, y,1e-8), 'Good market does not clear'
    assert np.isclose(24-h, l,0.0), 'Labor market does not clear'

    return p,h,y,pi,Inc,c,l

    #Print statement which are not used. 
    print(f'price = {p:.1f}')
    #print(f'consumer 1 consumption of good 1 = {c11:.1f}')
    #print(f'consumer 1 consumption of good 2 = {c21:.1f}')
    #print(f'consumer 2 consumption of good 1 = {c12:.1f}')
    #print(f'consumer 2 consumption of good 2 = {c22:.1f}')

#Making an interactive plot of the solution
def interactive_figure(alpha2,e11,e21,e12,e22):
    
    alpha1_vec = np.linspace(0,1,10)
    
    p_vec = np.empty(len(alpha1_vec))
    c11_vec = np.empty(len(alpha1_vec))
    c21_vec = np.empty(len(alpha1_vec))
    c12_vec = np.empty(len(alpha1_vec))
    c22_vec = np.empty(len(alpha1_vec))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    # a. calculations
    for i, alpha1 in enumerate(alpha1_vec):
        p_vec[i],c11_vec[i],c21_vec[i],c12_vec[i],c22_vec[i] = consumption(alpha1,alpha2,e11,e21,e12,e22)
    
    # b. figure
    ax.plot(alpha1_vec, c11_vec, label='Consumption of good 1 by consumer 1')
    ax.plot(alpha1_vec, c21_vec, label='Consumption of good 2 by consumer 1')
    ax.plot(alpha1_vec, c12_vec, label='Consumption of good 1 by consumer 2')
    ax.plot(alpha1_vec, c22_vec, label='Consumption of good 2 by consumer 2')
    ax.set_xlim([0.1,0.9]) # 
    ax.set_ylim([0,20]) #
    ax.set_title("Distribution of goods")
    ax.legend(loc= 'upper right')

