# ENGSCI263: Forecasting Geothermal Impacts
# wellbore_model.py

# PURPOSE:
# Function library for the AUTOUGH2 wellbore model, its surrogate and plotting

#########################################################################################
#
# 	Functions:
#		forward_model: runs the forward model, which is either AUTOUGH2 or its surrogate
#		wellbore_model: runs the AUTOUGH2 model
#		surrogate_model: runs the surrogate model for the AUTOUGH2 model
#		plot_distributions: plots distributions of permebaility and porosity
# 		plot_ensemble: plots an ensemble of forward models
#
#########################################################################################

# import modules
import numpy as np
from t2grids import *
from t2data import *
from t2listing import *
# from iapws import IAPWS97
import csv
from bisect import bisect_left
import matplotlib.pyplot as plt

ts = 16. 		# text size for plotting

# an object for storing simulation output
class Output(object):
    def __init__(self):
        self.final_earnings = None
    def compute_earnings(self, price_electricity):
        """Compute total dollars earned by well over its life.

        Parameters:
        -----------
            price_electricity : float
                cost of 1 kWhr of electricity
        """
        # midpoints of time vector
        self.tmid = 0.5*(self.t[1:]+self.t[:-1])
        # step size of time vector
        dt = self.t[1:]-self.t[:-1]
        # midpoints of power generation vector
        Wgmid = 0.5*(self.Wg[:1]+self.Wg[:-1])
        # income per step
        income = Wgmid*dt*365.25*24*price_electricity
        # cumulateive income
        self.earnings = np.cumsum(income)/1.e6
        # final earnings
        self.final_earnings = self.earnings[-1]

        # compute pressure decline in the reservoir
        self.final_pressure_decline = 20.-(self.mw[-1]/self.mw[0]*(20-14.5)+14.5)
        
#########################################################################################
#	
# Do not modify this function.
#
#########################################################################################
def load_experiment():
    """Load experiment samples, normalize against range of surrogate model

    Returns:
    --------
        ks : array-like
            Vector of permeability samples.
        phis : array-like
            Vector of porosity samples.
    """
        # load in raw results from experiment
    ks, phis = np.genfromtxt('class_results.txt', delimiter = ',', skip_header = 1).T
    phis = phis*100.
        # range of surrogate model
    min_phi,max_phi = [1,10] 		# porosity
    mu_phi2,std_phi2 = 0.5*(min_phi+max_phi), 0.2*(max_phi-min_phi)
    min_k,max_k = [1.e-16, 5.e-16]	# permeability
    logmin_k,logmax_k = np.log10(min_k), np.log10(max_k)
    logmu_k2,logstd_k2 = 0.5*(logmin_k+logmax_k), 0.25*(logmax_k-logmin_k)
        # statistics of experimental results
    mu_phi,std_phi = np.mean(phis), np.std(phis)
    mu_k,std_k = np.mean(ks), np.std(ks)
    logmu_k,logstd_k = np.mean(np.log10(ks)), np.std(np.log10(ks))
        # transform experiment results to new range
    phis = (phis - mu_phi)/std_phi*std_phi2+mu_phi2
    ks = 10**((np.log10(ks) - logmu_k)/logstd_k*logstd_k2+logmu_k2)
    
    return ks, phis
        
#########################################################################################
#	
# Do not modify this function.
#
#########################################################################################
def forward_model(permeability, porosity, price_electricity, costs=0.7, plot = False, surrogate = False, show_me_more = False):
    """Run wellbore simulation and calculate the power output given 
        by a backpressure turbine connected to the well output

    Parameters:
    -----------
        permeability : float
            permeability (m^2)
        porosity : float
            porosity (%)
        price_electricity : float
            cost of 1 kWhr of electricity
        costs : float
            total cost of running well over ten year period
        plot : bool
            optional boolean to plot outcome of simulation
        surrogate : bool
            switch to run forward model using surrogate, otherwise AUTOUGH2 used
    
    Returns:
    --------
        out : Output
            model object with simulation parameters and outputs
    """
    # create Output object for storing simulation parameters and output
    out = Output()
    # save simulation parameters
    out.k = permeability
    out.phi = porosity
    # create time vector
    out.t = np.linspace(0., 10., 120+1)[1:]	# one value per month
    # choose forward model: surrogate or AUTOUGH2
    if surrogate: model = surrogate_model
    else: model = wellbore_model
    # run forward model
    out.Wg, out.mw, out.m2, out.hw, out.h1g, out.h2 = model(out.t, permeability, porosity)
    # calculate earnings from well over 10 year period	
    out.compute_earnings(price_electricity)
    # if no plot requested, return power output
    if not plot: return out
    
    # plotting
    fig = plt.figure(figsize = [10., 15.])	# open figure
    ax1 = plt.axes([0.10, 0.55, 0.8, 0.35])		
    ax1a = ax1.twinx()
    ax2 = plt.axes([0.10, 0.10, 0.8, 0.35])		
    ax2a = ax2.twinx()
    
    # first plot, reservoir takes
        # plotting
    ln1 = ax1.plot(out.t,out.mw, 'b-', label = 'q(t), mass take from reservoir')
    if show_me_more:
        ln2 = ax1.plot(out.t,out.m2, 'b--', label = 'steam mass at turbine')
        ln3 = ax1a.plot(out.t,out.hw, 'r-', label = 'enthalpy in reservoir')
    else:
        dP = out.mw/out.mw[0]*(20-14.5)+14.5
        ln1 = ax1a.plot(out.t,dP, 'b-')
        ax1a.axhline(19., color='b',linestyle=':')
        ax1a.text(5, 19.05, 'impact on surface features', va = 'bottom', size=ts)
        # legend
    #lns = ln1+ln2+ln3
    #lbs = [ln.get_label() for ln in lns]
    ax1.legend(loc=3,prop={'size':ts})
        # limits
    if show_me_more:
        ax1a.set_ylim([800,3000])	
    xlim = [0,10]
    ax1.set_xlim(xlim)
    ax1a.set_xlim(xlim)
        # labels
    if show_me_more:
        ax1a.plot(out.t,out.t*0+out.h2, 'r:')
        ax1a.plot(out.t,out.t*0+out.h1g, 'r:')
        ax1a.text(out.t[-1], out.h2, 'enthalpy leaving turbine ', size = ts, ha = 'right', va = 'bottom')
        ax1a.text(out.t[-1], out.h1g, 'enthalpy entering turbine ', size = ts, ha = 'right', va = 'bottom')
        ax1a.text(out.t[-1], out.hw[-1], 'enthalpy in the reservoir ', size = ts, ha = 'right', va = 'bottom')
        
    ax1.set_ylabel('flow rate (kg/s)', size = ts)
    if show_me_more:
        ax1a.set_ylabel('enthalpy (kJ/kg)', size = ts)
    else:
        ax1a.set_ylabel('reservoir pressure (MPa)', size = ts)
    
    # second plot, power output and income
        # plotting
    ln1 = ax2.plot(out.t, out.Wg, 'k-', label = 'power output')
    ln2 = ax2a.plot(out.tmid, out.earnings, 'k--',label = 'earnings')
        # legend
    lns = ln1+ln2
    lbs = [ln.get_label() for ln in lns]
    ax2.legend(lns,lbs,loc=8,prop={'size':ts})
        # limits
    ax2.set_xlim(xlim)
    ax2a.set_xlim(xlim)
        # labels
    ax2a.plot(out.t,out.t*0+costs, 'r:')
    ax2a.text(0.5*out.t[-1], costs, 'well construction costs', size = ts, ha = 'center', va = 'bottom')
    ax2.set_ylabel('power output (kW)', size = ts)
    ax2a.set_ylabel('earnings (million NZD)', size = ts)
    
    # plot upkeep
    axs = [ax1,ax1a,ax2,ax2a]
    for ax in axs:
        for t in ax.get_xticklabels()+ax.get_yticklabels(): t.set_fontsize(ts)
        ax.set_xlabel('time (yrs)', size = ts)
    
    # display figure
    plt.savefig('model.png')
    plt.show()
    return out
    
#########################################################################################
#	
# Do not modify this function.
#
#########################################################################################
def wellbore_model(t, permeability, porosity):
    """Run wellbore simulation with AUTOUGH2 and calculate the power output given 
    by a backpressure turbine connected to the well output

    Parameters:
    -----------
        t : array-like 
            time vector (yr)
        permeability : float
            permeability (m^2)
        porosity : float
            porosity (% )
    Returns:
    --------
        Wg : array-like
            Power output (kW)
        mw : array-like
            Mass extraction rate (kg/s)
        m2 : array-like
            Mass rate of steam into turbine (kg/s)
        hw : array-like
            Reservoir enthalpy (kJ/kg)
        h1g : array-like
            Enthalpy at turbine input (kJ/kg)
        h2 : array-like
            Enthalpy at turbine output (kJ/kg)
    """
    
    # Reservoir parameters
    p0 = 20.e6	# Reservoir pressure (Pa)
    T0 = 250.	# Reservoir temperature (C)
    h = 200.	# reservoir thickness (m)
    
    # Well parameters
    z = 1500.	# well length (m)
    rw = .5		# well radius (m)
    
    # Simulation parameters
    tmax = 10.*365*24*3600	# time simulated (10 years)(s)
    dt = 365*24*3600/12		# timestep (1 month) (s)
    
    # Turbine parameters
    p1 = 10.e5	# separation pressure (Pa)
    p2 = 1.e5	# turbine back pressure (Pa)
    effG = .97	# generator efficiency

    # build grid and rocktypes	
    list_r = [rw]					# list of radial thicknesses
    list_r.extend([10.]*10)
    list_r.extend(list(np.logspace(1, 2, num = 70)))
    grid = t2grid().radial(list_r, [h], convention=2)											# create radial grid
    grid.add_rocktype(rocktype('WELL ', permeability = [.3e-10]*3, porosity = .99))				# add well rocktype
    grid.add_rocktype(rocktype('SANDS', permeability = [permeability]*3, porosity = porosity/100.))	# add SANDS rocktype
    grid.delete_rocktype('dfalt')															# delete the dfalt rocktype
    grid.blocklist[0].rocktype = grid.rocktype['WELL ']									# assign WELL to the first block
    for blk in grid.blocklist[1:]:													# assign SANDS to every other block
        blk.rocktype = grid.rocktype['SANDS']

    # DAT file writing
    dat = t2data()			# open empty dat file
    dat.title = 'radial sandy reservoir'	# simulation name
    dat.simulator = 'AUTOUGH2.2EW'			# simulator
    dat.start = True	
    dat.grid = grid				# assign grid

    # time parameters
    dat.parameter['const_timestep'] = (t[1]-t[0])*365*24*3600	# assign timestep (s)
    dat.parameter['default_incons'] = [p0, T0]	# assign reservoir temperature and pressure
    dat.parameter['gravity'] = 0.				# assign gravity in reservoir (horizontal radial flow, not needed)
    dat.parameter['max_timesteps'] = 3500		# assign maximum number of timesteps
    dat.parameter['tstart'] = 0.				# assign start time (s)
    dat.parameter['tstop'] = t[-1]*365*24*3600	# assign end time (s)

    # generator
    PI = 2*np.pi*200.*permeability/(np.log(10./np.sqrt(np.pi)))		# Productivity index (m3)
    Pw = p1+(1.e3*9.81*z)*.9		# equivalent separator pressure at reservoir depth,.
    # = separation pressure + weight of the water column*coeff
    # coeff: part of the water will flash to steam, so the wellbore isn't completely filled with liquid water
    # a coupled wellbore simulator would have to be implemented to get a more accurate estimation
    dat.add_generator(t2generator(		# Add well
        name = 'WEL 1', 				# name of the well
        block = ' a  1', 				# well position
        type = 'DELG',					# type: here on deliverability 
        ltab = 1, 						# number of layers producing
        gx = PI,						# set productivity index
        ex = Pw,					# set wellbore pressure
        ))

    # save and run
    dat.write('TEST.DAT')		# write AUTOUGH2 .DAT file
    dat.run(simulator='AUTOUGH2_5.exe', silent = True)	# run AUTOUGH2

    # get results
    lst = t2listing('TEST.LISTING', skip_tables = ['connection'])	# extract listing file
    [(th, h), (tm, m)] = lst.history([('g', (' a  1', 'WEL 1'), 'Enthalpy'), ('g', (' a  1', 'WEL 1'), 'Generation rate')])
    th = np.array(th)*(1/24.)*(1/3600.)*(1/365.)	# convert time from seconds to years
    h = np.array(h)*1.e-3				# convert enthalpy from J/kg to kJ/kg
    m = np.array(m)*-1					# switch to positive flow rates
    
    # extract history values in order to get a unified time vector
    mw = np.array([])
    hw = np.array([])
    for time in t:		# loop through time output vector
        index = (np.abs(th-time)).argmin()	# get the index of the closest time from the history vector
        mw = np.append(mw, m[index])	# get the mass flow rate at the index
        hw = np.append(hw, h[index])	# get the enthalpy at the flow rate
    
    # fixed thermodynamics variables (only for p1 = 10.e5 Pa and p2 = 1.e5 Pa )
    h1l = 762.7		# saturated liquid specific enthalpy at 10 bar (kJ/kg)
    h1g = 2777.1	# saturated steam specific enthalpy at 10 bar (kJ/kg)
    h2l = 417.4		# saturated liquid sspecific enthalpy at 1 bar (kJ/kg)
    h2g = 2674.9	# saturated steam specific enthalpy at 1 bar (kJ/kg)

    s1g = 6.585		# saturated steam specific entropy at 10 bars (kJ/kg.K)
    s2l = 1.303		# saturated liquid specific entropy at 1 bar (kJ/kg.K)
    s2g = 7.359		# saturated steam specific entropy at 1 bar (kJ/kg.K)

    # thermodynamics calculations (depend only on turbine parameters)
    x2is = (s1g - s2l)/(s2g - s2l)	# dryness at turbine output considering an isentropic transformation
    h2is = x2is*(h2g-h2l)+h2l		# enthalpy at turbine output considering an isentropic transformation
    At = .425*(h1g-h2is)			# Baumann rule coefficient (kJ/kg)
    h2 = (h1g-At*(1-h2l/(h2g-h2l)))/(1+At/(h2g-h2l))	# h at turbine output considering the real transformation (kJ/kg)
    
    # thermodynamics calculations (depend on well output)
    x1 = (hw - h1l)/(h1g - h1l)	# dryness at p1
    m2 = mw * x1				# turbine mass flow rate (kg.m-3)
    Wg = m2*(h1g-h2)*effG		# Power output of the generator (kW)
    
    return Wg, mw, m2, hw, h1g, h2

#########################################################################################
#	
# Do not modify this function.
#
#########################################################################################	
def surrogate_model(t, permeability, porosity):
    """Run surrogate model for wellbore_model()
    
     Parameters:
     -----------
        t : array-like
            time vector (yr)
        permeability : float
            permeability (m^2)
        porosity : float
            porosity (%)
    Returns:
    --------
        Wg : array-like
            Power output (kW)
        mw : array-like
            Mass extraction rate (kg/s)
        m2 : array-like
            Not populated, for consistency only
        hw : array-like
            Not populated, for consistency only
        h1g : array-like
            Not populated, for consistency only
        h2 : array-like
            Not populated, for consistency only
    """

    # read csv files containing calculated parameters of the surrogate model for different permeability porosity combinations
    list_key = ['a', 'b', 'c', 'd']	# vector listing parameters name
    table = {}						# empty dictionary
    for key in list_key:			# for each parameter
        table[key] = np.loadtxt(open(key+'.csv','rb'), delimiter = ',')	# extract parameters values from csv file
    k = np.loadtxt(open('list_k.csv','rb'), delimiter = ',')		# vector of permeability values for which the surrogate model parameters have been calculated
    phi = np.loadtxt(open('list_phi.csv','rb'), delimiter = ',')	# vector of porosity values for which the surrogate model parameters have been calculated
    # bilinear interpolation
    j = bisect_left(k, permeability)	# index of closest higher permeability for which the surrogate model parameters have been calculated
    i = bisect_left(phi, porosity)		# index of closest higher porosity for which the surrogate model parameters have been calculated
    dk = np.array([k[j]-permeability, permeability-k[j-1]])	# vector for bilinear interpolation calculation
    
    try:
        dphi = np.array([phi[i]-porosity, porosity-phi[i-1]])	# vector for bilinear interpolation calculation
    except:
        print(phi, porosity, i)
    #dphi *= 100.
    theta = []	# empty vector
    for key in list_key:	# for each parameter
        square_matrix = np.transpose(table[key][np.ix_([i-1,i],[j-1,j])])	# matrix containing the 4 parameter values framing the chosen permeability porosity combination
        theta.append(dk.dot(square_matrix.dot(dphi))/((k[j]-k[j-1])*(phi[j]-phi[j-1]))) # calculation of the parameter value for the chosen permeability porosity combination

    # surrogate model
    a,b,c,d = theta
    h1l = 762.7		# saturated liquid specific enthalpy at 10 bar (kJ/kg)
    h1g = 2777.1	# saturated steam specific enthalpy at 10 bar (kJ/kg)
    h2l = 417.4		# saturated liquid sspecific enthalpy at 1 bar (kJ/kg)
    h2g = 2674.9	# saturated steam specific enthalpy at 1 bar (kJ/kg)

    s1g = 6.585		# saturated steam specific entropy at 10 bars (kJ/kg.K)
    s2l = 1.303		# saturated liquid specific entropy at 1 bar (kJ/kg.K)
    s2g = 7.359		# saturated steam specific entropy at 1 bar (kJ/kg.K)

    # thermodynamics calculations (depend only on turbine parameters)
    x2is = (s1g - s2l)/(s2g - s2l)	# dryness at turbine output considering an isentropic transformation
    h2is = x2is*(h2g-h2l)+h2l		# enthalpy at turbine output considering an isentropic transformation
    At = .425*(h1g-h2is)			# Baumann rule coefficient (kJ/kg)
    h2 = (h1g-At*(1-h2l/(h2g-h2l)))/(1+At/(h2g-h2l))	# h at turbine output considering the real transformation (kJ/kg)
    effG = .97	# generator efficiency

    Wg = np.array([a*np.exp(-b*time)-c*time+d for time in t])

    m2 = Wg/(effG*(h1g-h2))
    x1 = (1100 - h1l)/(h1g - h1l)
    mw = m2/x1

    return [Wg, mw, None, None, None, None]

#########################################################################################
#	
# Do not modify this function.
#
#########################################################################################
def plot_distributions(porosity, permeability, poro_stats, perm_stats):
    """Plot fitted distributions of porosity and permeability

    Args:
        porosity (numpy.array): sampled porosity values
        permeablity (numpy.array): sampled permeability values
        poro_stats (list): two item list comprising mean and standard deviation of fitted porosity distribution
        perm_stats (list): two item list comprising mean and standard deviation of fitted permeability distribution
    """
    # plotting
    fig = plt.figure(figsize = [10., 7.])	# open figure
    ax1 = plt.axes([0.10, 0.10, 0.38, 0.80])		
    ax2 = plt.axes([0.55, 0.10, 0.38, 0.80])		
    
    # first plot, porosity distribution
        # calculations
    h, e = np.histogram(porosity, bins = 10)
    h = h/(np.sum(h)*1.*(e[1]-e[0]))		# normalize to PDF
    phi = np.linspace(0,10,101) 				# plot vector
        # only plot normal distribution if parameters supplied
    if None in poro_stats:
        Pphi = 0.*phi
    else:
        mu,std = poro_stats						# unpack stats vector
        Pphi = np.exp(-(phi-mu)**2/(2*std**2))/np.sqrt(2*std**2*np.pi)	# compute PDF
        # plotting
    ax1.bar(e[:-1], h, e[1]-e[0], color = [0.5, 0.5, 0.5])
    ax1.plot(phi,Pphi, 'k-', lw = 2)
        # limits
    ax1.set_xlim([0,10])
        # labels
    ax1.set_xlabel('porosity', size = ts)
    ax1.set_ylabel('PDF', size = ts)
    
    # second plot, permeability distribution
        # calculations
    permeability = np.log10(permeability)
    h, e = np.histogram(permeability, bins = 10)
    h = h/(np.sum(h)*1.*(e[1]-e[0]))
        # only plot normal distribution if parameters supplied
    if None not in perm_stats:
        mu,std = perm_stats
        k = np.linspace(mu-3*std,mu+3*std,101)				
        Pk = np.exp(-(k-mu)**2/(2*std**2))/np.sqrt(2*std**2*np.pi) # compute PDF
        ax2.plot(k, Pk, 'k-', lw = 2)
        # plotting
    ax2.bar(e[:-1], h, e[1:]-e[:-1], color = [0.5, 0.5, 0.5])
        # labels
    ax2.set_xlabel('log10(permeability)', size = ts)
    ax2.set_ylabel('PDF', size = ts)
    ax2.set_xticks([-16,-15.8,-15.6,-15.4])
    ax2.set_xlim([-16.1, -15.25])
    #ax2.set_xscale('log')
    
    # plot upkeep
    for ax in [ax1,ax2]:
        for t in ax.get_xticklabels()+ax.get_yticklabels(): t.set_fontsize(ts)
        ax.set_ylabel('PDF', size = ts)
    
    plt.savefig('distributions.png')
    plt.show()
    

#########################################################################################
#	
# Do not modify this function.
#
#########################################################################################
def plot_ensemble(outs, costs):
    """Plot ensemble of models

    Parameters:
    -----------
        outs : array-like
            list of model Output objects to plot results for
        costs : float
            cost of drilling and operating well for 10 years (million NZD)
    """
    N = len(outs)
    
    # plotting initiation and parameters
    print('\nPlotting...')
    N_bins = 100	# number of bins for hist plots
    alpha_plot = 1./np.sqrt(N)	# opacity of plots
    fig = plt.figure(figsize = [8., 12])	# open figure
    
    ax1 = plt.axes([0.10,0.72,0.35,0.25])
    ax2 = plt.axes([0.55,0.72,0.35,0.25])
    ax5 = plt.axes([0.32,0.40,0.35,0.25])
    #ax6 = plt.axes([0.55,0.40,0.35,0.25])
    ax3 = plt.axes([0.10,0.08,0.35,0.25])
    ax4 = plt.axes([0.55,0.08,0.35,0.25])
    
    # plot k distribution
    nk = ax2.hist([out.k for out in outs], bins = int(np.sqrt(N)), facecolor='green', alpha=0.75)[0]	# hist plot + highest number of samples in one bin
    xlim = ax2.get_xlim(); ax2.set_xlim(xlim)
    ax2.text(.8*(xlim[1]-xlim[0])+xlim[0], max(nk), str(N)+' samples', fontsize = ts, ha = 'right', va = 'bottom')	# show sample number
    
    # plotting upkeep
    ax2.set_ylim(0., 1.2*max(nk))
    ax2.set_xlabel('Permeability ('+r'$m^2$)', fontsize = ts)
    ax2.tick_params(labelsize=ts)
    
    # phi distribution
    nphi = ax1.hist([out.phi for out in outs], int(np.sqrt(N)), facecolor='green', alpha=0.75)[0]	# hist plot + highest number of samples in one bin
    xlim = ax1.get_xlim(); ax1.set_xlim(xlim)
    ax1.text(.8*(xlim[1]-xlim[0])+xlim[0], max(nphi), str(N)+' samples', fontsize = ts, ha = 'right', va = 'bottom')	# show sample number
    
    # plotting upkeep
    ax1.set_ylim(0., 1.2*max(nphi))
    ax1.set_xlabel('Porosity (%)', fontsize = ts)
    ax1.tick_params(labelsize=ts)
    
    # power output
    for out in outs:	# for each sample
        ax3.plot(out.t, out.Wg, c='k', alpha = alpha_plot, linewidth = 1.)	# show power output
        ax5.plot(out.t, out.mw/out.mw[0]*(20.-14.5)+14.5, c='b', alpha = alpha_plot, linewidth = 1.)	# show pressure decline
    
    # plotting upkeep
    ax3.set_xlim(out.t[0], out.t[-1])
    ax3.set_ylim(0., 80.)
    ax3.set_xlabel('Time (years)', fontsize = ts)
    ax3.set_ylabel('Power output (kW)', fontsize = ts)
    ax3.tick_params(labelsize=ts)
    
    # cumulative energy produced / money earned
    for out in outs:	# for each sample
        ax4.plot(out.tmid, out.earnings, c='k', alpha = alpha_plot, linewidth = 1.)	# show cumulated money earned
    ax4.plot([0., 10.], [costs, costs], 'r--')	# show profitability threshold
    ax4.text(.3, costs*1.03, 'Profitability threshold', fontsize = ts)
    ax5.text(.3, 18.95, 'Impact threshold', fontsize = ts, va = 'top')
    ax5.axhline(19., color='r',linestyle='--')
    
    # plotting upkeep
    ylabel = 'Money earned (million NZD)'
    ax4.set_xlim(0.,10.)
    ax4.set_ylim([0., 0.8])
    ax4.set_xlabel('Time (years)', fontsize = ts)
    ax4.set_ylabel(ylabel, fontsize = ts)
    ax4.tick_params(labelsize=ts)
    
    ax5.set_xlim(out.t[0], out.t[-1])
    ax5.set_xlabel('Time (years)', fontsize = ts)
    ax5.set_ylabel('Reservoir pressure (MPa)', fontsize = ts)
    ax5.tick_params(labelsize=ts)
    
    plt.savefig('ensemble.png')	
    plt.show()
    