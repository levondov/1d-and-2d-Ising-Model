from math import *
import matplotlib.pyplot as plt
import matplotlib.pylab as pyl
import numpy as np
import time
from scipy.stats import norm
from scipy import optimize
from scipy.optimize import leastsq
import matplotlib as mpl
import matplotlib.animation as manimation

cmap2 = mpl.colors.ListedColormap(['C0','C1'])


def delta_E(zx,zy,chainz):
    """
    2d ising model (Probability Density Function)
    """
    chain_N = len(chainz)

    if zx < (chain_N-1) and zy < (chain_N-1):
    # no edge problem
        H = chainz[zx,zy]*chainz[zx+1,zy] + chainz[zx,zy]*chainz[zx-1,zy] + chainz[zx,zy]*chainz[zx,zy+1] + chainz[zx,zy]*chainz[zx,zy-1] + h*chainz[zx,zy]
    else:
        if zx == 0 and zy == 0:
            H = chainz[zx,zy]*chainz[zx+1,zy] + chainz[zx,zy]*chainz[chain_N-1,zy] + chainz[zx,zy]*chainz[zx,zy+1] + chainz[zx,zy]*chainz[zx,chain_N-1] + h*chainz[zx,zy]
        elif zx == (chain_N-1) and zy == (chain_N-1):
            H = chainz[zx,zy]*chainz[0,zy] + chainz[zx,zy]*chainz[zx-1,zy] + chainz[zx,zy]*chainz[zx,0] + chainz[zx,zy]*chainz[zx,zy-1] + h*chainz[zx,zy]
        elif zx == 0 and zy == (chain_N-1):
            H = chainz[zx,zy]*chainz[zx+1,zy] + chainz[zx,zy]*chainz[chain_N-1,zy] + chainz[zx,zy]*chainz[zx,0] + chainz[zx,zy]*chainz[zx,zy-1] + h*chainz[zx,zy]
        elif zx == (chain_N-1) and zy == 0:
            H = chainz[zx,zy]*chainz[0,zy] + chainz[zx,zy]*chainz[zx-1,zy] + chainz[zx,zy]*chainz[zx,zy+1] + chainz[zx,zy]*chainz[zx,chain_N-1] + h*chainz[zx,zy]
        else:
            if zx == 0:
                H = chainz[zx,zy]*chainz[zx+1,zy] + chainz[zx,zy]*chainz[chain_N-1,zy] + chainz[zx,zy]*chainz[zx,zy+1] + chainz[zx,zy]*chainz[zx,zy-1] + h*chainz[zx,zy]
            elif zx == (chain_N-1):
                H = chainz[zx,zy]*chainz[0,zy] + chainz[zx,zy]*chainz[zx-1,zy] + chainz[zx,zy]*chainz[zx,zy+1] + chainz[zx,zy]*chainz[zx,zy-1] + h*chainz[zx,zy]
            elif zy == 0:
                H = chainz[zx,zy]*chainz[zx+1,zy] + chainz[zx,zy]*chainz[zx-1,zy] + chainz[zx,zy]*chainz[zx,zy+1] + chainz[zx,zy]*chainz[zx,chain_N-1] + h*chainz[zx,zy]
            elif zy == (chain_N-1):
                H = chainz[zx,zy]*chainz[zx+1,zy] + chainz[zx,zy]*chainz[zx-1,zy] + chainz[zx,zy]*chainz[zx,0] + chainz[zx,zy]*chainz[zx,zy-1] + h*chainz[zx,zy]
            else:
                print 'this should never happen'
    return H

    if z == (chain_N-1):
        H = chainz[z-1]*chainz[z] + chainz[z]*chainz[0] + h*chainz[z]
    elif z == 0:
        H = chainz[chain_N-1]*chainz[z] + chainz[z]*chainz[z+1] + h*chainz[z]
    else: # no edge effect
        H = chainz[z-1]*chainz[z] + chainz[z]*chainz[z+1] + h*chainz[z]
    return H

def periodic_step(xstep,K):
    """
    Take care of periodic boundary conditions
    """
    if xstep >= n[K]:
        # longer than chain
        return (xstep-n[K])
    elif xstep < 0:
        # shorter than chain length
        return (xstep+n[K])
    else:
        # not a periodic boundary condition
        return xstep

def measure_E(chainz,K):
    """
    Measure Energy of the system
    """
    temp=0
    for i in range(1,n[K]-2):
        for j in range(1,n[K]-2):
            temp += chainz[i,j]*chainz[i+1,j] + chainz[i,j]*chainz[i-1,j] + chainz[i,j]*chainz[i,j+1] + chainz[i,j]*chainz[i,j-1] + h*chainz[i,j]
    return -1*temp  
    
def measure_m(chainz):
    return sum(sum(chainz[1:-1,1:-1]))
           
fig = plt.figure("main plot",figsize=(9,10))

###################### Parameters ##################################
# trials
N = 1

# beta,T,H
kb = 1.0

T = np.repeat(3.5,N) # or np.linspace(0.1,5,N)
beta = 1./(kb*T)
hs = np.repeat(0,N) # np.linspace(-1.5,1.5,N)

# iterations of steps
L = int(1e5)
# chain length
n = [int(50)]*N # nxn
# random step length
steps = np.repeat(n[0]/2,N) # or np.linspace(1,10,N)

# Thermalization time (need to adjust manually by looking at results
T_time = 50

# chain initial condition, 1 = all spins down, 2 = all spins random (up or down)
chain_start = 1

# Measurables
measure_rate = 50 # rate to measure parameters in Metropolos time (i.e. every x Metropolos steps)
bmeasure_E = False # boolean for measuring Energy every measure_rate
bmeasure_m = True # boolean for measuring magnetization every measure_rate
bmeasure_E_avg = False # boolean for measuring average Energy every iteration
bmeasure_m_avg = False # boolean for measuring average magnetization every iteration

# extra plots and calculations
plot_grid = True
if plot_grid:
    plt.ion()
    fig_grid = plt.figure('grid plot')


# plot average magnetization vs magnetic field
plot_mag_vs_field = False # make sure bmeasure_m_avg is turned on

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test')
writer = FFMpegWriter(fps=120, metadata=metadata)
      
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
    chain = np.random.randint(chain_start,size=(n[j],n[j])) #random chain of 1d spins [-1,0)U(0,1]
    chain = chain + (chain-1) # turn list of 0&1 to -1&1
        
    # initial starting point random
    x = np.random.randint(0,n[j])
    y = np.random.randint(0,n[j])
        
    # define measurables
    vx = [] # chain idx selected by random jumps
    vy = [] # chain idy selected by random jumps
    m = [] # magnetization
    E = [] # Energy
    vx.append(x)
    vy.append(y)
    m.append(measure_m(chain))
    E.append(measure_E(chain,j))
    
    with writer.saving(plt.figure('grid plot'),"2dising.mp4",100):
        for l in range(L): # number of iterations        
            # proposal
            x_can = x + np.random.randint(-step,step)
            y_can = y + np.random.randint(-step,step)
            x_can = periodic_step(x_can,j) # check for periodic effects
            y_can = periodic_step(y_can,j)
            vx.append(x_can)
            vy.append(y_can)
            chain_can = np.copy(chain)
            chain_can[x_can,y_can] = -1*chain_can[x_can,y_can] 
            # acceptance prob
            prob = min([1.,exp(beta[j]*delta_E(x_can,y_can,chain_can))/exp(beta[j]*delta_E(x_can,y_can,chain))])
            u = pyl.uniform(0,1)
            if u < prob:
                acceptance_rate += 1
                x = x_can
                y = y_can
                chain = chain_can
                
            # measure variables
            if l % measure_rate == 0: # measure every x number of iterations
                if bmeasure_m:
                    m.append(measure_m(chain))
                if bmeasure_E:
                    E.append(measure_E(chain,j))

            if plot_grid and l % measure_rate == 0:
                print 'Metropolos step: ' + str(l) + ' | Magnetization: ' + str(m[-1])
                plt.figure('grid plot')
                # tell imshow about color map so that only set colors are used
                plt.cla()
                img = plt.imshow(chain,cmap=cmap2)
                cb = plt.colorbar(img,ticks=[-1,1],cmap=cmap2,boundaries=[-1,0,1])
                writer.grab_frame()
                plt.pause(0.05)
                cb.remove()

    # magnetization average after each run
    if bmeasure_m_avg:
        m_avg.append(np.mean(m[T_time:])/n[j]**2)
    # Energy average after each run
    if bmeasure_E_avg:
        E_avg.append(np.mean(E[T_time:])/n[j]**2)

    # plotting results
    plt.figure("main plot")
    plt.subplot(3,2,1)
    plt.cla()
    plt.hist(vx,bins=n[j],edgecolor='black',label='Chain length = '+str(n[j]))
    plt.xlabel('x_proposed chain loc for each Metropolos step')
    #plt.title('2D Ising Model run with h=0, J/Kb=1, T=3.5, and acceptance= ~0.5')
    plt.legend()
    plt.subplot(3,2,2)
    plt.cla()
    plt.hist(vy,bins=n[j],edgecolor='black',label='Chain length = '+str(n[j]))
    plt.xlabel('y_proposed chain loc for each Metropolos step')
    plt.legend()
    plt.subplot(3,2,(3,4))
    plt.cla()
    plt.axvline(T_time,color='r',label='Thermalization')
    plt.plot(m,label='_nolegend_')
    plt.xlabel('')
    plt.ylabel('Net Magnetization')
    plt.legend()
    plt.grid(True)
    plt.subplot(3,2,(5,6))
    plt.cla()
    plt.plot(E)
    plt.axvline(T_time,color='r',label='Thermalization')
    plt.legend()
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
    s = np.linspace(25,1000,1000)

    # Plotting stuff
    plt.figure()
    plt.scatter(M_time,error,marker='.',color='C0',label='Simulation')
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
    #hss = np.linspace(hs[0],hs[-1],200)
    T = T[0] # temperature is an array so take first value
    #m_theo = ( np.sinh(hss/T) + (np.sinh(hss/T)*np.cosh(hss/T)/np.sqrt(
    #    np.sinh(hss/T)**2 + np.exp(-4/T))) )/ ( np.cosh(hss/T) + 
    #    np.sqrt(np.sinh(hss/T)**2 + np.exp(-4/T)) )
    plt.figure()
    plt.scatter(hs,m_avg_sizes[0],marker='.',label='Sim L = 10x10')
    plt.scatter(hs,m_avg_sizes[1],marker='.',label='Sim L = 25x25')
    plt.scatter(hs,m_avg_sizes[2],marker='.',label='Sim L = 80x80')         
    #plt.plot(hss,m_theo,color='C1',label='Theoretical')
    plt.xlabel('Magnetic field H')
    plt.ylabel('Average Magnetization')
    plt.title('Average Magnetization vs Magnetic Field')
    plt.legend()

# Magnetization vs Temperature
if 0:
    plt.figure()
    plt.scatter(T,m_avg)
    plt.xlabel('Temperature')
    plt.ylabel('Magnetization')

plt.show()        
