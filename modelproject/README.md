# Model analysis project

Our project is titled **Koopman Economy with CES utility function**. We model a production economy. First, we find the analytical solution in our baseline model, where the utility function is a Cobb-Douglas function. Secondly, we find the solution numerical by an algorithm described in the notebook. The extension of our model is to use a CES utility function in order to allow different degrees of substitutability, as the Cobb-Douglas function has strong assumption about income and substitution effect. 

The **results** of the project can be seen from running [modelproject.ipynb](modelproject.ipynb). We find that higher degree of substitutability will lead to lower price effects of change in the supply (modelled through a change in $\beta$). This is because the consumer can easier substitute toward leisure, if the goods are closer to be perfect complements. 

**Dependencies:** Apart from a standard Anaconda Python 3 installation and packages installed through out the course, the project requires no further installations.