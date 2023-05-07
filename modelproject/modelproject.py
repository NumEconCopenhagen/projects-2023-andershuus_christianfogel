from scipy import optimize
import numpy as np
import ipywidgets as widgets # Interactive plots
import matplotlib.pyplot as plt


#Analytical solution, i.e. section 2
def consumption(alpha1,alpha2,e11,e21,e12,e22): 
    """
        
    Args:

        alpha1: relative preferences for consumer 1
        alpha2: relative preferences for consumer 2
        e11: Initial endowment of good 1 for consumer 1
        e21: Initial endowment of good 2 for consumer 1
        e12: Initial endowment of good 1 for consumer 2
        e22: Initial endowment of good 2 for consumer 2
        
    Returns:
    
        c11 = Consumption of good 1 for consumer 1
        c21 = Consumption of good 2 for consumer 1    
        c12 = Consumption of good 1 for consumer 2    
        c22 = Consumption of good 2 for consumer 2
        p = relative price of good 1  
        
    """
    # Initally find the price
    p = (alpha1*e21+alpha2*e22)/(e11*(1-alpha1)+e12*(1-alpha2))

    #Then define the income for each consumer: 
    I_1=p*e11+e21
    I_2=p*e12+e22

    #Consumption of good 1 and good 2 for consumer 1
    c11 = alpha1*I_1/p
    c21 = (1-alpha1)*I_1

    #Consumption of good 1 and good 2 for consumer 2
    c12 = alpha2*I_2/p
    c22 = (1-alpha2)*I_2

    #Check that supply is equal to demand for both goods. Demand cannot be smaller due to assumption about monotonicity
    assert np.isclose(e11+e12, c11+c12,0.0), 'supply is not equal to demand for good 1'
    assert np.isclose(e21+e22, c21+c22,0.0), 'supply is not equal to demand for good 2'

    return p,c11,c21,c12,c22

    #Print statement which are not used. 
    #print(f'price = {p:.1f}')
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

