import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from math import exp, log, sqrt, pi, fsum
from ipywidgets import interact, FloatSlider, Dropdown

class Markov_Widget():
    
    def __init__ (self):
        interact(self.solve_and_plot,
                 P1 = FloatSlider(value=45, min=0, max=100, step=1, description='P1', continuous_update=False),
                 P2 = FloatSlider(value=20, min=0, max=100, step=1, description='P2', continuous_update=False),
                 P3 = FloatSlider(value=65, min=0, max=100, step=1, description='P3', continuous_update=False),
                 P4 = FloatSlider(value=50, min=0, max=100, step=1, description='P4', continuous_update=False),
                 P5 = FloatSlider(value=20, min=0, max=100, step=1, description='P5', continuous_update=False),
                 P6 = FloatSlider(value=15, min=0, max=50, step=1, description='P6', continuous_update=False),
                 P7 = FloatSlider(value=1.0, min=0, max=5, step=0.1, description='P7', continuous_update=False),
                 P8 = FloatSlider(value=0.023, min=0, max=2, step=0.01, description='P8', continuous_update=False),
                 P9 = FloatSlider(value=29, min=0, max=100, step=1, description='P9', continuous_update=False),
                 P10 = FloatSlider(value=15, min=0, max=50, step=1, description='P10', continuous_update=False),
                 P11 = FloatSlider(value=2e-5, min=0, max=1, step=1e-6, description='k1', continuous_update=False),
                 P12 = FloatSlider(value=2e-5, min=0, max=1, step=1e-6, description='k2', continuous_update=False),
                 P13 = FloatSlider(value=0.5, min=0, max=3, step=0.1, description='gK', continuous_update=False))

    def solve_and_plot(self,P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13):
        
        # The parameter vector
        
        SSA_data = np.loadtxt("SS.txt",dtype='float')
        V = SSA_data[:,0]
        I = SSA_data[:,1]        
        init_params = [P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13]
        step_length = 1000
        
        dats = Activation(init_params,V,step_length)
        model_I = dats['I_peak']
        Po = dats['Po']
        t = dats['t']
                
        plt.plot(V,I,'b-')
        plt.plot(V,model_I,'r-')
        plt.xlabel('Step potential [mV]')
        plt.ylabel('Current (A/F)')
        plt.legend((r'Experiment',r'Model'))
        plt.show()        
        
def f(t,y,P):
# A Markov state model for the rapidly activating inward rectifying potassium current (IKur)      
    
    V = P[-1] # For simplicity we pass voltage as a parameter (when would this be a bad idea?)
    dim = y.shape
    
    alpha = exp((V-P[0])/P[1])
    beta = exp((V-P[2])/P[3])*exp(-(V+P[4])/P[5])/(P[6]+P[7]*exp(-(V+P[8])/P[9]))
    dy = np.zeros((dim))
    #states
    O = int(0)
    I = int(1)
    C4 = int(2)
    C3 = int(3)
    C2 = int(4)
    C1 = int(5)
    B = int(6)
    
    n = 7 # number of states
    A = np.zeros((n,n)) # initialize rate mass matrix for output

    #calculate transitions between closed states
    A[C1, C2] = beta
    A[C2, C1] = 4*alpha
    A[C2, C3] = 2*beta
    A[C3, C2] = 3*alpha
    A[C3, C4]  = 3*beta
    A[C4, C3]  = 2*alpha
    
    #calculate transitions between closed and open states
    A[O, C4]  = alpha
    A[C4, O]  = 4*beta
    
    #calculate transitions between inactive and open states
    A[O, I]  = P[10] # k1
    A[I, O]  = P[11] # k2
    
    #calculate transitions between open and blocked states
    A[O,B] = 0
    A[B,O] = 0
    
    for i in range(n):
        A[i,i] = -fsum(A[:,i])
    
    dy = A.dot(y)
    return [dy]

def Activation(P,data,duration):
    
    #configure inputs
    V = data
    ntsteps = int(duration*10)
    
    #initialize outputs
    t = np.zeros((ntsteps))
    y = np.zeros((ntsteps,7))
    model_peaks = np.zeros((len(V)))
    Po_out = np.zeros((ntsteps,len(V)))
    
    #step through the test potentials in your reference data
    for n,i in enumerate(V):
        
        # First reset your initial conditions to the holding potential at the beginning of each step
        # Question: Is this always a good idea? Why or why not?
        
        V_H = -70
        y0 = Init(V_H, P)
        t0 = 0
        
        # collect the parameters and include the step voltage in the last position
        if n==0 and len(P)<14:        
            P.append(i)
        else:
            P[-1]=i
            
        # integrate the model at the current step potential                         
        r = ode(f).set_integrator('vode', method='bdf', with_jacobian=False, rtol=1e-5, max_step=0.1)
        r.set_initial_value(y0,t0).set_f_params(P)
        dt = 0.1
        ind = 0
        while r.successful() and r.t < duration and ind < ntsteps:
            r.integrate(r.t+dt)
            t[ind]=r.t
            y[ind,:]=r.y
            ind = ind+1
        
        #store the open probability data and calculate the error metrics
        Po_out[:,n] = y[:,0]
        model_peaks[n] = P[-2]*y[:,0].max()
        out = {'t':t, 'Po':Po_out,'I_peak':model_peaks}        
        
    return out

def Init(V,P):
    
    y0 = np.zeros((7))
    alpha = exp((V-P[0])/P[1])
    beta = exp((V-P[2])/P[3])*exp(-(V+P[4])/P[5])/(P[6]+P[7]*exp(-(V+P[8])/P[9]))
    gamma = beta/alpha
    Kblk = P[10]/P[11]
    o = 1/((1+gamma)**4+Kblk)
    
 #    The initial conditions for the model:
    y0 = [o, o*Kblk, o*4*gamma, o*6*gamma**2, o*4*gamma**3, o*gamma**4, 0]
    
    return y0


