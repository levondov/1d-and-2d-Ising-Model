from math import *
import matplotlib.pyplot as plt
import matplotlib.pylab as pyl
import numpy as np
import time
from scipy.stats import norm
from scipy import optimize
from scipy.optimize import leastsq

def delta_E(z,chainz):
    """
    1d ising model (Probability Density Function)
    """
    chain_N = len(chainz)
    
    if z == (chain_N-1):
        H = chainz[z-1]*chainz[z] + chainz[z]*chainz[0] + h*chainz[z]
    elif z == 0:
        H = chainz[chain_N-1]*chainz[z] + chainz[z]*chainz[z+1] + h*chainz[z]
    else: # no edge effect
        H = chainz[z-1]*chainz[z] + chainz[z]*chainz[z+1] + h*chainz[z]
    return H

def periodic_step(xstep):
    """
    Take care of periodic boundary conditions
    """
    if xstep >= n:
        # longer than chain
        return (xstep-n)
    elif xstep < 0:
        # shorter than chain length
        return (xstep+n)
    else:
        # not a periodic boundary condition
        return xstep

def measure_E(chainz):
    """
    Measure Energy of the system
    """
    temp=0
    it = iter(chainz)
    for k in it:
        temp = temp + k*next(it) + h*k
    return -1*temp  
    
def measure_m(chainz):
    return sum(chainz[1:-1])    
           
fig = plt.figure(figsize=(8,9))

###################### Parameters ##################################
# trials
N = 50

# beta,T,H
kb = 1.0

T = np.repeat(1.9,N) # or np.linspace(0.1,5,N)
beta = 1./(kb*T)
hs = np.repeat(0,N) # or np.linspace(-1.5,1.5,N)

# iterations of steps
L = int(1e5)
# chain length
n = int(1e3) 
# random step length
steps = np.repeat(10,N) # or np.linspace(1,10,N)

# Thermalization time (need to adjust manually by looking at results
T_time = 200

# chain initial condition, 1 = all spins down, 2 = all spins random (up or down)
chain_start = 1

# Measurables
measure_rate = 100 # rate to measure parameters in Metropolos time (i.e. every x Metropolos steps)
bmeasure_E = True # boolean for measuring Energy every measure_rate
bmeasure_m = True # boolean for measuring magnetization every measure_rate
bmeasure_E_avg = False # boolean for measuring average Energy every iteration
bmeasure_m_avg = False # boolean for measuring average magnetization every iteration

# extra plots and calculations

# plot average magnetization vs magnetic field
plot_mag_vs_field = False # make sure bmeasure_m_avg is turned on
      
####################################################################
m_avg = []
E_avg = []
        
for j in range(N): # number of trials
    # setup initial parameters
    acceptance_rate = 0 # measures acceptance rate
    run_time = time.time() # measures run time on computer
    step = steps[j]
    h = hs[j]
    
    # create chain
    chain = np.random.randint(chain_start,size=n) #random chain of 1d spins [-1,0)U(0,1]
    chain = chain + (chain-1) # turn list of 0&1 to -1&1
        
    # initial starting point random
    x = np.random.randint(0,n)
        
    # define measurables
    v = [] # chain idx selected by random jumps
    m = [] # magnetization
    E = [] # Energy
    v.append(x)
    m.append(measure_m(chain))
    E.append(measure_E(chain))
    
    for l in range(L): # number of iterations        
        # proposal
        x_can = x + np.random.randint(-step,step)
        x_can = periodic_step(x_can) # check for periodic effects
        v.append(x_can)
        chain_can = np.copy(chain)
        chain_can[x_can] = -1*chain_can[x_can]
            
        # acceptance prob
        prob = min([1.,exp(beta[j]*delta_E(x_can,chain_can))/exp(beta[j]*delta_E(x_can,chain))])
        u = pyl.uniform(0,1)
        if u < prob:
            acceptance_rate += 1
            x = x_can
            chain = chain_can
            
        # measure variables
        if l % measure_rate == 0: # measure every x number of iterations
            if bmeasure_m:
                m.append(measure_m(chain))
            if bmeasure_E:
                E.append(measure_E(chain))

    # magnetization average after each run
    if bmeasure_m_avg:
        m_avg.append(np.mean(m[T_time:])/n)
    # Energy average after each run
    if bmeasure_E_avg:
        E_avg.append(np.mean(E[T_time:])/n)
    
    # plotting results
    plt.subplot(3,1,1)
    plt.cla()
    plt.hist(v,bins=n/10,edgecolor='black',label='Chain length = '+str(n))
    plt.xlabel('proposed chain location for each Metropolos step')
    plt.title('1D Ising Model run with h=0, J/Kb=1, T=1.9, and acceptance= ~0.5')
    plt.legend()
    plt.subplot(3,1,2)
    plt.cla()
    plt.axvline(T_time,color='r',label='Thermalization')
    plt.plot(m,label='_nolegend_')
    plt.xlabel('')
    plt.ylabel('Net Magnetization')
    plt.legend()
    plt.grid(True)
    plt.subplot(3,1,3)
    plt.cla()
    plt.plot(E)
    plt.xlabel('every ' + str(measure_rate) + ' Metropolos steps')
    plt.ylabel('Energy')
    plt.grid(True)
    plt.pause(0.05)            
        
        
        
    print 'trial ' + str(j) + ' | acceptance: ' + str(acceptance_rate/float(L)) + ' | proposal size: ' + str(step) + ' | Run time: ' + str(time.time()-run_time)[0:5]

##### Extra calculations and plotting here #####
# Need to slightly edit them manually
#
# Energy and error calculations
if 0:
    # energy calculations
    E_ave = []
    error = []
    M_time = []

    # calculate error and Energy measurement
    for i in range(T_time+25,T_time+1000):
        E_ave.append(np.mean(E[T_time:i]))
        error.append(np.std(E[T_time:i])/(i-T_time))
        M_time.append(i-T_time+10)

    # fit function 1/sqrt(N)
    # Fit the first set
    fitfunc = lambda p, x: p[0]*1/np.sqrt(x) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    p0 = [1] # Initial guess for the parameters
    p1, success = optimize.leastsq(errfunc, p0[:], args=(M_time, error))
    s = np.linspace(15,1000,1000)

    # Plotting stuff
    plt.figure()   
    plt.scatter(M_time,error,marker='.',label='Experiment')
    plt.plot(s,fitfunc(p1,s),color='C1',label='Theoretical ($y=\\frac{1}{\sqrt{N}}$)')
    plt.xlabel('Metropolos time after Thermalization (every 100 steps)')
    plt.ylabel('Error')
    plt.legend()
    plt.title('Error in Energy Measurement after Thermalization')
    plt.figure()
    plt.plot(M_time,E_ave)
    plt.xlabel('Metropolos time after Thermalization (every 100 steps)')
    plt.ylabel('Energy')
    plt.title('Measure of Energy vs Metropolos time')
    
# Magnetization vs Field calculations
if plot_mag_vs_field:
    hss = np.linspace(hs[0],hs[-1],200)
    T = T[0] # temperature is an array so take first value
    m_theo = ( np.sinh(hss/T) + (np.sinh(hss/T)*np.cosh(hss/T)/np.sqrt(
        np.sinh(hss/T)**2 + np.exp(-4/T))) )/ ( np.cosh(hss/T) + 
        np.sqrt(np.sinh(hss/T)**2 + np.exp(-4/T)) )
    plt.figure()
    plt.scatter(hs,m_avg,marker='.',label='Simulation')
    plt.plot(hss,m_theo,color='C1',label='Theoretical')
    plt.xlabel('Magnetic field H')
    plt.ylabel('Average Magnetization')
    plt.title('Average Magnetization vs Magnetic Field for L = '+str(n))
    plt.legend()

# Magnetization vs Temperature
if 0:
    plt.figure()
    plt.scatter(T,m_avg)
    plt.xlabel('Temperature')
    plt.ylabel('Magnetization')

plt.show()        
