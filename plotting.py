
#########################################################################################
#
# Function library for the AUTOUGH2 wellbore model
#
# 	Functions:
#		plot_sampling: runs wellbore model (or its surrogate) for a given number of permeability porosity samples
#		and display the results
#
#########################################################################################

# import modules and functions
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import datetime
from wellbore_model import *

def plot_sampling(N, sample_k, sample_phi, price_electricity, pthres):
    """Plot wellbore model results for N samples

    Args:
        N (integer): number of samples
        sample_k (callable): permeability distribution function
        sample_phi (callable): porosity distribution function
        price_electricity (float): electricity cost (NZD/kWh)
        pthres (float): profitability threshold (NZD)
    """
    
    # output number of samples chosen on screen
    print str(N)+' model simulations'
    
    # model running parameters
    t = np.linspace(0., 10., 120+1)[1:]	# one value per month
    N_switch = 5	# iteration from which the surrogate model is used
    
    # save grid parameters in csv, and load them
    list_k = np.loadtxt(open('list_k.csv','rb'), delimiter = ',')	# vector of permeability values for which the surrogate model parameters have been calculated
    list_phi = np.loadtxt(open('list_phi.csv','rb'), delimiter = ',')	# vector of porosity values for which the surrogate model parameters have been calculated
    range_k = [list_k[0], list_k[-1]]	# possible range of permeability values
    range_phi = [list_phi[0], list_phi[-1]]	# possible range of porosity values
    
    # simulations
    samples = []	# initiate save list of dictionaries for samples
    i = 0			# initiate iteration number
    dt_wellbore = 0.	# initiate total runtime wellbore model 
    dt_surrogate = 0.	# initiate total runtime surrogate model 
    while i < N:	# check if sample number is reached
        k = sample_k(range_k)			# get a permeability sample
        phi = sample_phi(range_phi)		# get a porosity sample
        if range_k[0] < k < range_k[1] and range_phi[0] < phi < range_phi[1]:	# check if permeability and porosity are in range
            samples.append({})		# add a dictionary to the save list
            samples[i]['k'] = k		#save permeability
            samples[i]['phi'] = phi	# save porosity
            if i < N_switch:	# check if surrogate model should be run
                t0 = datetime.datetime.now()		# time beginning simulation
                samples[i]['w'] = wellbore_model(t, k, phi)	# run wellbore model
                t1 = datetime.datetime.now()		# time end simulation
                dt = (t1-t0).total_seconds()		# time difference
                dt_wellbore += dt					# add simulation runtime to total runtime wellbore model 
                print 'Wellbore model number '+str(i+1)+' done'
            else:
                if i == N_switch:	# check if this is the first surrogate model to be used
                    print 'Switch to surrogate model'
                t0 = datetime.datetime.now()		# time beginning simulation
                samples[i]['w'] = surrogate_model(t, k, phi)	# run surrogate model
                t1 = datetime.datetime.now()		# time end simulation
                dt = (t1-t0).total_seconds()		# time difference
                dt_surrogate += dt					# add simulation runtime to total runtime surrogate model 
            samples[i]['g'] = np.cumsum(samples[i]['w'])*365*2*1.e-3*price_electricity	# calculate cumulated money earned from energy generation
            i += 1	# update sample number

    # display time difference wellbore/surrogate model
    print '\nSimulations done'
    if N > N_switch:	# check if surrogate model has been used
        s1 = dt_wellbore+dt_surrogate	# label1
        s2 = seconds=dt_wellbore*N/5	# label2
        label1 = '\nTotal run time (first '+str(N_switch)+' samples use wellbore model, next samples use surrogate model): '
        label2 = '\nExtrapolated run time if surrogate was not used: '
        for s, label in zip([s1, s2], [label1, label2]):	# better time display
            hours, remainder = divmod(s, 3600)
            minutes, seconds = divmod(remainder, 60)
            s_print = ''
            if hours > 0:
                s_print += str(int(hours))+' hours, '
            if minutes > 0 or hours > 0:
                s_print += str(int(minutes))+' minutes, '+str(int(seconds))+' seconds'
            else:
                s_print += str(seconds)+' seconds'
            print label
            print s_print
    
    # plotting initiation and parameters
    print '\nPlotting...'
    N_bins = 100	# number of bins for hist plots
    text_size = 16.	# text size in plot
    alpha_plot = 1./np.sqrt(N)	# opacity of plots
    fig = plt.figure(figsize = [20., 13.])	# open figure
    
    # plot k distribution
    ax1 = plt.subplot2grid((4,4), (0,0), colspan=2)		# open subplot
    bins = np.linspace(range_k[0], range_k[1], N_bins)	# create bins
    nk = ax1.hist([s['k'] for s in samples], bins, facecolor='green', alpha=0.75)[0]	# hist plot + highest number of samples in one bin
    ax1.text(.8*(range_k[1]-range_k[0])+range_k[0], max(nk), str(N)+' samples', fontsize = text_size)	# show sample number
    
    # plotting upkeep
    ax1.set_xlim(range_k)
    ax1.set_ylim(0., 1.2*max(nk))
    ax1.set_xlabel('Permeability ('+r'$m^2$)', fontsize = text_size)
    ax1.tick_params(labelsize=text_size)
    
    # phi distribution
    ax1 = plt.subplot2grid((4,4), (1,0), colspan=2)			# open subplot
    bins = np.linspace(range_phi[0], range_phi[1], N_bins)	# create bins
    nphi = ax1.hist([s['phi'] for s in samples], bins, facecolor='green', alpha=0.75)[0]	# hist plot + highest number of samples in one bin
    ax1.text(.8*(range_phi[1]-range_phi[0])+range_phi[0], max(nphi), str(N)+' samples', fontsize = text_size)	# show sample number
    
    # plotting upkeep
    ax1.set_xlim(range_phi)
    ax1.set_ylim(0., 1.2*max(nphi))
    ax1.set_xlabel('Porosity (%)', fontsize = text_size)
    ax1.tick_params(labelsize=text_size)
    
    # power output
    ax1 = plt.subplot2grid((4,4), (0,2), colspan=2, rowspan=2)	# open subplot
    for s in samples:	# for each sample
        ax1.plot(t, s['w'], c='k', alpha = alpha_plot, linewidth = 1.)	# show power output
    
    # plotting upkeep
    ax1.set_xlim(t[0], t[-1])
    ax1.set_ylim(0., 80.)
    ax1.set_xlabel('Time (years)', fontsize = text_size)
    ax1.set_ylabel('Power output (kW)', fontsize = text_size)
    ax1.tick_params(labelsize=text_size)
    
    # cumulated energy produced / money earned
    ax1 = plt.subplot2grid((4,4), (2,0), colspan=2, rowspan=2)	# open subplot
    for s in samples:	# for each sample
        ax1.plot(t, s['g'], c='k', alpha = alpha_plot, linewidth = 1.)	# show cumulated money earned
    ax1.plot([0., 10.], [pthres, pthres], 'r--')	# show profitability threshold
    ax1.text(.3, pthres*1.03, 'Profitability threshold', fontsize = text_size)
    
    # plotting upkeep
    ylim = 800.
    ylabel = 'Money earned (kNZD)'
    ax1.set_xlim(0.,10.)
    ax1.set_ylim(0., ylim)
    ax1.set_xlabel('Time (years)', fontsize = text_size)
    ax1.set_ylabel(ylabel, fontsize = text_size)
    ax1.tick_params(labelsize=text_size)
    
    # distribution energy produced
    ax1 = plt.subplot2grid((4,4), (2,2), colspan=2, rowspan=2)	# open subplot
    list_G = [s['g'][-1] for s in samples]	# list of money earned after 10 years for each samples
    list_G1 = [G for G in list_G if G < pthres]		# list of money earned after 10 years below profitability
    list_G2 = [G for G in list_G if G >= pthres]	# list of money earned after 10 years above profitability
    range_G = [0., 800.]	# range of bins
    bins = np.linspace(range_G[0], range_G[1], int((range_G[1]-range_G[0])/10.)+1)		# create bins
    if len(list_G1) > 0 :	# check if samples with money earned after 10 years below profitability exist
        n1  = max(plt.hist(list_G1, bins, facecolor='red', alpha=0.75)[0])	# maximum number of samples in one bin for those below profitability
    else: n1 = 0.
    if len(list_G2)>0:	# check if samples with money earned after 10 years above profitability exist
        n2 = max(plt.hist(list_G2, bins, facecolor='green', alpha=0.75)[0])	# maximum number of samples in one bin for those above profitability
    else: n2 = 0.
    ylim = 1.2*max(n1, n2) # y upper boundary for plotting purpose
    ax1.plot([pthres, pthres], [0., ylim], 'k--')	# show profitability threshold
    pctage_G1 = 100.*len(list_G1)/len(list_G)	# percentage of samples profitable
    pctage_G2 = 100.*len(list_G2)/len(list_G)	# percentage of samples non profitable
    ax1.text(pthres/2, .9*ylim, 'No profit after 10 years\n('+str(pctage_G1)+'%)',fontsize=text_size,ha='center')
    ax1.text((pthres+range_G[1])/2,.9*ylim,'Profit after 10 years\n('+str(pctage_G2)+'%)',fontsize=text_size,ha='center')
    
    # add a 'best fit' line using a normal distribution
    if len(list_G)>9:	# check if enough samples simulated
        mu, sigma = np.mean(list_G), np.std(list_G)	# get mean and standard deviation
        if sigma > 0.:	# check if standard deviation above 0
            y = mlab.normpdf(bins, mu, sigma)*(bins[1]-bins[0])*len(list_G)	# calculate normal distribution 
            l = ax1.plot(			# show best fit
                bins, y, 'b--', linewidth=1, 
                label = r'$\mu$: '+str(round(mu, 2))+' kNZD\n'+r'$\sigma$: '+str(round(sigma, 2))+' kNZD'
                )
            ax1.legend(loc = 'center right', fontsize = text_size, framealpha = 0.)
    
    # plotting upkeep
    ax1.set_ylim(0., ylim)
    ax1.set_xlabel('Money earned after 10 years (kNZD)', fontsize = text_size)
    ax1.tick_params(labelsize=text_size)
    
    # save and show
    plt.tight_layout()
    plt.savefig('lab3_plot1_'+str(N)+'.png', bbox_inches = 'tight')			



    
                                                                                                                                                                                            