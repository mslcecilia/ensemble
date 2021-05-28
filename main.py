# ENGSCI263: Forecasting Geothermal Impacts
# main.py

# PURPOSE:
# To USE a prior distribution to CONSTRUCT an ensemble of models.

# PREPARATION:
# Notebook uncertainty.ipynb.

# SUBMISSION:
# Show the instructor that you can replicate Fig. 5 in the lab document.

# INSTRUCTIONS:
# Jump to section "if __name__ == "__main__":" at bottom of this file.

# import modules and functions
from time import time
import numpy as np
from matplotlib import pyplot as plt
from wellbore_model import *


#########################################################################################
#	
# Task 1: Model familiarization
#
#########################################################################################
def get_familiar_with_model():
    ''' This function runs the wellbore model ONCE for given parameters.
    '''
    # **to do**
    # 1. uncomment and choose values for k and phi below
    # 2. research the price of electricity (NZD per kWhr)
    # 3. run FORWARD_MODEL and inspect the output
    # 4. ANSWER the questions in the lab document
    
    # reservoir properties 
    #k = 			# permeability (m^2)
    #phi =  		# volume (%) 
    
    # economic parameters
    #price_electricity =  # NZD per kWhr 
    costs = 0.70 			# costs (million NZD) of drilling and running a well for 10 years
    
    # run the forward model and plot the outcome
    out = forward_model(k, phi, price_electricity, costs, plot = True, show_me_more=False)
    print("Total earnings over ten years is %3.2f million NZD"%out.final_earnings)
    print("Reservoir pressure drop after ten years is %2.1f MPa"%out.final_pressure_decline)

    # WHAT IS 'out'?
    # A model output object. 
    # Its attributes contain information about the model run.
    # e.g., out.k stores the permeability
    #       out.final_earnings stores the total income from the well after ten years
    #       out.final_pressure_decline stores the total pressure drop in the reservoir
    #          after ten years

    # EXTENSION
    # For those especially interested...
    # in 'forward_model', set show_me_more=True and rerun this task
    # i.   Hot liquid water in the reservoir is converted to a mixture of steam (blue dashed
    #      line in top plot) and hot water at the surface. Why?
    # ii.  Enthalpy in the reservoir is lower than enthalpy at the surface. Why?
    # iii. Enthalpy drops as steam passes through a turbine. Why?
    return

#########################################################################################
#	
# Task 2: Construct prior distributions
#
#########################################################################################
def construct_priors():
    ''' This function plots prior distributions over data.
    '''
    # **to do**
    # 1. uncomment and CHOOSE values for mean and standard deviations of phi and logk
    # 2. ANSWER the questions in the lab document

    # import results from gingernut experiment + calibration assignment
    ks, phis = load_experiment()

    # parameters for porosity distribution
    mu_phi = None			
    std_phi = None			

    # parameters for permeability distribution
    mu_logk = None			
    std_logk = None 		

    # plot permeability and porosity samples as histograms, overlay fitted distributions
    plot_distributions(phis, ks, [mu_phi, std_phi], [mu_logk, std_logk])	
    return

#########################################################################################
#	
# Task 3: Open fun_with_randoms.py and complete the exercises.
#
#########################################################################################

#########################################################################################
#	
# Task 4: Create an ensemble of models
#
#########################################################################################
def model_ensemble(N):
    ''' This function runs the wellbore model N times.

        Parameters:
        -----------
        N : int
            Number of times to run the model.

        Returns:
        --------
        outs : array-like
            List of model output objects.
    '''
    # **to do**
    # 1. copy-paste your prior parameters from TASK 2
    # 2. CREATE random samples of phi and k (see TASK 3)
    # 3. RUN forward_model (see TASK 1) and save the output you are interested in
    # 4. ANSWER the questions in the lab document
    
    # load in the data
    ks, phis = load_experiment()

    # limits of surrogate model
    min_phi,max_phi = [1,10] 		# porosity
    min_k,max_k = [1.e-16, 5.e-16]	# permeability

    # 1. copy-paste your prior parameters from above
    # mu_phi = 
    # std_phi = 
    # mu_logk = 
    # std_logk = 
    
    # economic parameters
    price_electricity = 0.25 # NZD per kWhr 

    n = 5 			# number of models to run "properly" before switching to the (fast) surrogate
    outs = []		# empty list for saving model output objects
    t0 = time() 	# start time
    cnt = 0
    while cnt < N:
        # 2a. use function np.random.randn() to create a sample	porosity
        # phi = 
        
        # 2b. use function np.random.randn() to create a sample	permeability
        # note, because we have used a normal prior for LOG of k, sampling
        # is a two step process: (i) generate a random number logk, (ii) convert to k
        # logk = 
        # k = 
        
        # if values outside range of surrogate model, reject and restart loop
        if ((k<min_k) or (k>max_k)): 
            continue
        if ((phi<min_phi) or (phi>max_phi)): 
            continue
            
            # use surrogate or full model? (see lab document: surrogate models)
        if cnt > n-1: 
            surrogate = True
        else: 
            surrogate = False
        
        # print model iteration 
        if not surrogate:
            print("Running model %i"%(cnt+1))
            
        # print timing information
        if cnt == n:
            t1 = time() # time at which switch to surrogate occurs
            print("Switching to surrogate model")
            m, s = divmod(t1-t0, 60); h, m = divmod(m, 60)
            print("%i models run with simulator took %d:%02d:%02d"%(n,h, m, s))
            m, s = divmod(N*(t1-t0)/n, 60); h, m = divmod(m, 60)
            print("Extrapolated time to run %i models with simulator is %d:%02d:%02d"%(N, h, m, s))
            
        # 3. run the forward model and save earnings and pressure decline
        # (refer to Task 1, lines 49-51, set plot=False, and surrogate=surrogate)
        # out = 
        # outs.append(
        
        # iterate counter
        cnt +=1
        
    # print timing information
    t2 = time() 	# time when all simulations finished
    m, s = divmod(t2-t0, 60); h, m = divmod(m, 60)
    print("Actual time using surrogate is %d:%02d:%02d"%(h, m, s))
        
    # plot results of model ensemble
    plot_ensemble(outs, costs=0.7)

    return outs

#########################################################################################
#	
# Task 5: Quantify forecast outcomes
#
#########################################################################################
def decision(predictions, threshold, xlabel='predictions', thresholdlabel='threshold'):
    ''' This function plots a histogram of the model predictions and tests these
        against a threshold value.

        Parameters:
        -----------
        predictions : array-like
            List of model predictions (e.g., final earnings, reservoir pressure decline).
        threshold : float
            A particular prediction value, where the goal is to stay below or above.
        xlabel : string
            Quantity for which a prediction has been made. Used as label for x-axis.
        thresholdlabel : string
            Quantity for which threshold corresponds. Used as label for vert. line.

        Notes:
        ------
        The plotting functionality is complete. You should complete the calculation that
        computes the probability to exceed the threshold.
    '''
    # **to do**
    # 1. WRITE code to compute the probability of exceeding a threshold prediction
    # 2. ANSWER the questions in the lab document

    # create figure and axes
    fig = plt.figure(figsize = [10., 7.5])	# open figure
    ax = plt.axes([0.15,0.15,0.70,0.70])
    
    # plot a histogram of predictions 
    N = len(predictions)
    h,e = np.histogram(predictions, int(np.sqrt(N)))
    ax.bar(e[:-1], h, e[1]-e[0], label='model ensemble')

    # set axis limits and add threshold
    ax.axvline(threshold, color='r', linestyle=':', label=thresholdlabel)
    
    # labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel('frequency')

    ax.legend()
    
    # display plot
    plt.show()

    # 1. compute likelihood to exceed threshold
    # *hint* use a for loop, or look up 'np.where'
    # Pexceed =  
    
    print("Probability to exceed threshold %3.2f"%Pexceed)

if __name__ == "__main__":
    # Comment/uncomment each of the functions below as you complete the tasks
    
    # TASK 1: Read the instructions in the function definition.
    get_familiar_with_model()
    
    # TASK 2: Read the instructions in the function definition.
    #construct_priors()

    # TASK 3: Open the file fun_with_randoms.py and complete the tasks

    # TASK 4: Read the instructions in the function definition.
    # leave the command below UNCOMMENTED when completing TASKS 5 and 6
    #N = 100
    #outs = model_ensemble(N)
    
    # TASK 5: Read the instructions in the function definition
    # leave the TASK 4 command UNCOMMENTED when completing TASK 5
    #earnings = [out.final_earnings for out in outs]
    #costs = 0.7
    #decision(predictions=earnings, threshold=costs, xlabel='final earnings', thresholdlabel='well costs')

    # TASK 6: Copy-paste-modify the commands from TASK 5 above to test
    #    whether pressure decline is less than the 1 MPa threshold.
    # leave the TASK 4 command UNCOMMENTED when completing TASK 6
    #declines = [
    #decision(
    