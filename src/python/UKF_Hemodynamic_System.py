import numpy as np
from scipy.linalg import cholesky
import math
from numpy.random import randn
import matplotlib.pyplot as plt
#from filterpy.stats import plot_covariance_ellipse
from scipy.integrate import solve_ivp
import matplotlib.patches as mpatches
import sys


class SigmaWeights:

    def __init__(self, x, P, alpha, beta, kappa):
        self.x = x
        self.P = P
        self.n = len(x)
        self.alpha = alpha
        self.beta = beta
        self.kappa = kappa

        self.sigmaPoints = np.zeros((2*self.n+1,self.n))
        self.weights = np.zeros(2*self.n+1)

        self._lambda = (alpha**2)*(self.n+kappa) - self.n
    

    def computeSigmaPoints(self):
        self.U = cholesky((self.n+self._lambda)*self.P)
        self.sigmaPoints[0] = self.x
        for i in range(self.n):
            self.sigmaPoints[i+1] = self.x + self.U[i]
            self.sigmaPoints[self.n+i+1] = self.x - self.U[i]
    

    def computeWeights(self):
        self.Wc = np.full(2*self.n+1, 1./(2*(self.n + self._lambda)))
        self.Wm = np.full(2*self.n+1, 1./(2*(self.n + self._lambda)))
        self.Wc[0] = self._lambda/(self.n+self._lambda) + (1. - alpha**2+beta)
        self.Wm[0] = self._lambda / (self.n +self._lambda)
    

    
def stateTransition(sigmas, f, time, dt):

    sigmasState = np.zeros_like(sigmas.sigmaPoints)
    for i, s in enumerate(sigmas.sigmaPoints):
        Rt=np.log2((16000*c)/(Rcref+Rdref))
        p_old=s[0]
        Rc=s[1]
        C=s[2]
        Rd=s[3]
        Zero=s[4]
           
        
        sol= solve_ivp(lambda t,y:f(t, y, C ,Rd ,Rc),
                       (time, time+dt), np.array([p_old]), method='LSODA')
        p=sol.y[0]
        p_new = p[-1]
        sigmasState[i,0] = p_new
        sigmasState[i,1] = Rc
        sigmasState[i,2] = C
        sigmasState[i,3] = Rd
        sigmasState[i,4] = Rc+Rd-Rt
        

    return sigmasState

def unscentedTransform(sigmas, sigmasState, Q):

    mean = np.dot(sigmas.Wm, sigmasState)
    covariance = np.zeros((sigmas.n,sigmas.n))
    for i in range(2*sigmas.n+1):
        temp = (sigmasState[i].reshape((sigmas.n,1))-mean.reshape((sigmas.n,1)))
        temp2 = np.matmul( temp, np.transpose(temp) )
        covariance += sigmas.Wc[i]*temp2

    covariance += Q
    


    return mean, covariance

    
def h_func(x):

    return np.array([x[0],x[4]])


def computeMeasurementSigmas(h, sigmasState):
    
    sigmasMeasurement = []
    for i in range(sigmasState.shape[0]):
        sigmasMeasurement.append(h(sigmasState[i,:]))

    return np.array(sigmasMeasurement)

    
def computeMeasurementMeanAndCovariance(sigmas, sigmasZ, R):

    lenColZ = sigmasZ.shape[1]
    mean = np.dot(sigmas.Wm, sigmasZ)
    covariance = np.zeros((lenColZ,lenColZ))
    for i in range(2*sigmas.n+1):
        temp = (sigmasZ[i].reshape((lenColZ,1))-mean.reshape((lenColZ,1)))
        temp2 = np.matmul( temp, np.transpose(temp) )
        covariance += sigmas.Wc[i]*temp2

    covariance += R
    
    return mean, covariance


def computeCrossCovariance_Pxz(sigmas, meanPredict, meanZ, sigmasState, sigmasZ):

    lenColZ = sigmasZ.shape[1]
    covariance = np.zeros((sigmas.n,lenColZ))
    for i in range(2*sigmas.n+1):
        temp1 = (sigmasState[i].reshape((sigmas.n,1))-meanPredict.reshape((sigmas.n,1)))
        temp2 = (sigmasZ[i].reshape((lenColZ,1))-meanZ.reshape((lenColZ,1)))
        temp3 = np.matmul( temp1, np.transpose(temp2) )
        covariance += sigmas.Wc[i]*temp3

    return covariance


def computeKalmanGain(crossCov, covZ):

    return np.dot(crossCov,np.linalg.inv(covZ))


def mkNoisyMeasurements(f, y0, time, R, std_constrain):
    
    C = np.log2(((2.5e-5)/c)/Cref)
    Rd = np.log2((13000*c)/Rdref)
    Rc = np.log2((1600*c)/Rcref)

    sol = solve_ivp(lambda t,y:f(t,y,C,Rd,Rc),
                    (time[0],time[-1]), np.array([y0]), t_eval=time, method='LSODA')
    p = sol.y[0]
    zP = p+np.random.normal(0,R,len(p))
    zR = np.random.normal(0,std_constrain,len(p))
    z = np.array([zP,zR])

    
    return p, z


def computeFinalSol(f, y0, time, *args):
    
    C = args[1]
    Rd = args[2]
    Rc = args[0]

    sol = solve_ivp(lambda t, y:f(t,y,C,Rd,Rc), (time[0],time[-1]), np.array([y0]), t_eval=time, method='LSODA')
    p = sol.y[0]
    
    return p


def f(t, y, C ,Rd ,Rc):
    
    t = t%0.8
    
    return 1/C*(Q(t)+Rc*C*dQ(t)-y/Rd+(Q(t))*Rc/Rd)

def Q(t):
    t=t%0.8
    S = 0.065
    return (1.4/(S*np.sqrt(np.pi)))*np.exp(-((t-0.36)**2)/(2*S**2))+5

def dQ(t):
    t=t%0.8
    dt = 1e-8
    return (Q(t+dt)-Q(t))/dt



if __name__ == "__main__":
    
    global c,Cref,Rcref,Rdref,Pref
    
    c = 0.000750061683  #unit conversion 1 dyn/cm^2=0.0007750061 mmHg
    Pref=50
    Cref=(2e-5)/c
    Rcref=1500*c
    Rdref=12000*c
    Rt=np.log2((500*c)/Rcref)+np.log2((25000*c)/Rdref)

    x = np.array([np.log2(80/Pref),np.log2((500*c)/Rcref) ,np.log2((0.00065/c)/Cref),np.log2((25000*c)/Rdref),np.log2((500*c)/Rcref)+np.log2((25000*c)/Rdref)])

    
    P = np.diag([8,8,8,8,8])

    
    
    alpha = 0.2
    beta = 2.
    kappa = 0
    R_adjust = 0.005
    z_std = 0.1
    std_constrain = 0.001
    #0.005
    cycles = 10
    np.random.seed(3)
    dt = 1/1000 #2*np.pi/100.
    time = np.arange(0, 0.8*cycles+dt, dt)
    #y0 = np.log2(78/Pref)
    y0=1
    sol, z  = mkNoisyMeasurements(f, y0, time, z_std, std_constrain)
    std = np.std(z)
    
    R = np.diag([R_adjust, 0.0000850])

    Qp = dt**2
    Qprocess =np.array([[ Qp, Qp, Qp, Qp, Qp],
                        [ Qp, Qp, Qp, Qp, Qp],
                        [ Qp, Qp, Qp, Qp, Qp],
                        [ Qp, Qp, Qp, Qp, Qp],
                        [ Qp, Qp, Qp, Qp, Qp]])*600


    '''
    Qprocess =np.array([[(dt**2),      0.0,      0.0,      0.0,      0.0],
                        [    0.0,  (dt**2),      0.0,      0.0,      0.0],
                        [    0.0,      0.0,  (dt**2),      0.0,      0.0],
                        [    0.0,      0.0,      0.0,  (dt**2),      0.0],
                        [    0.0,      0.0,      0.0,      0.0,  (dt**2)]])*650    
    '''
    predictions = []
    state = []
    res = []
    Filter_Order=[]
#    state.append(x) Uncomment if you want to show initial guess. Size of arrays will shift by 1 
    
    
    

    for i in range(len(time)):
        try:
            print('time = {}'.format(time[i]))
            Filter_Order.append(P[0,0])
#            Pplot=np.array([[P[0,0],P[0,1]],[P[1,0],P[1,1]]])
#            plot_covariance_ellipse((x[0],x[1]),Pplot, fc='g', alpha=0.2,)
#           plt.gca().grid(b=False);       
            sigmas = SigmaWeights(x, P, alpha, beta, kappa)
            
            
            sigmas.computeSigmaPoints()
            sigmas.computeWeights()       
            sigma_point_locations=np.array([[sigmas.sigmaPoints[:,0]],[sigmas.sigmaPoints[:,1]]])
#            plt.plot(sigma_point_locations[0,:],sigma_point_locations[1,:],'o')
            sigmasState = stateTransition(sigmas,f, time[i], dt)
            meanPredict, covPredict = unscentedTransform(sigmas, sigmasState, Qprocess)
            predictions.append(meanPredict)
            sigmasZ = computeMeasurementSigmas(h_func, sigmasState)
            meanZ, covZ = computeMeasurementMeanAndCovariance(sigmas, sigmasZ, R)
            crossCov = computeCrossCovariance_Pxz(sigmas, meanPredict, meanZ, sigmasState, sigmasZ)
            K = computeKalmanGain(crossCov, covZ)
            residual = np.subtract(z[:,i], meanZ)
            res.append(residual)
            x = np.add( meanPredict, np.dot(K,residual) )
            P = np.subtract( covPredict, np.dot(K, np.dot(covZ,np.transpose(K))) )
            state.append(x)
        except:
            try:
                test = linalg.inv(covPredict)
            except:
                print("Q is NOT invertable!!!")
                print(covPredict)
                sys.exit()
        
        
    plt.show()
    
    Fil_Ord=np.sqrt(np.array(Filter_Order))
    predictions = np.array(predictions)
    res = np.array(res)
    state = np.array(state)
    Rc = np.mean(Rcref/c*(2**(state[800:-1,1])))
    C = np.mean(Cref*c*(2**(state[800:-1,2])))
    Rd = np.mean(Rdref/c*(2**(state[800:-1,3])))
    state_normal=Pref*2**(state)
    print(Rc)
    print(C)
    print(Rd)
    
    finalSol = computeFinalSol(f, 80, time, Rc*c, C/c, Rd*c)
    realSol = computeFinalSol(f, 80, time, 1600*c, 2.5e-5/c, 13000*c)
    

#------------------Plot of pressure Vs time in Transformed Space--------------
    
    
    plt.figure(1, figsize=(10,13))
    
    
    plt.subplot(2,1,1)

    #plt.plot(time, sol, 'r')
    #green_patch = mpatches.Patch(color='r', label='Process Model')
    
    #plt.plot(time, z[0],'o', 'b')
    #blue_patch = mpatches.Patch(color='b', label='Pressure Measurements')
    
    plt.plot(time, state[:,0], 'g')
    red_patch = mpatches.Patch(color='g', label='Filter Estimate')
    
    #plt.legend(handles=[green_patch,blue_patch,red_patch])
    plt.legend(handles=[red_patch])
    
    plt.title('Pressure vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Transformed Pressure')

#---------------------Plot of Pressure vs time---------------------------
    

    
    plt.subplot(2,1,2)
    
    plt.plot(time, realSol, 'black')
    black_patch = mpatches.Patch(color='black', label='ODE solved with True Parameters')
    
    
    #plt.plot(time, finalSol, 'm', label='ODE sol with UKF parameters')
    #magenta_patch = mpatches.Patch(color='m', label='ODE sol with UKF parameters')
    

    #plt.legend(handles=[black_patch,magenta_patch])
    plt.legend(handles=[black_patch])
    
    plt.title('Pressure vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (mmHg)')
    
    plt.show()


#----------------------------Capacitance and Resistance------------------------
    
    
    
    
 #   plt.figure(2, figsize=(10,13))

    plt.figure(2, figsize=(16,12))
    
    plt.subplot(2,2,1)
    plt.plot(time, Rcref/c*(2**(state[:,1])),'g')
    plt.plot(time,1600*np.ones_like(time),'--')
 
    plt.title('Resistance Rc vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Resistance (dyn-s/cm^3)')
    
#--------------------------------------------------------------------------   
    plt.subplot(2,2,2)
    plt.plot(time, Cref*c*(2**(state[:,2])),'r')
    plt.plot(time, 2.5e-5*np.ones_like(time),'--')
    plt.title('C')
    
    plt.title('Capacitance Rd vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Capacitance (cm^3/dyn)')
    
    
    plt.subplot(2,2,3)
    plt.plot(time, Rdref/c*(2**(state[:,3])),'b')
    plt.plot(time, 13000*np.ones_like(time),'--')
    
    
    plt.title('Resistance Rd vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Resistance (dyn-s/cm^3)')
    
    plt.show()
    
#---------------------------Residuals---------------------------------------

    
    fig, ax1 = plt.subplots(1, 1, sharex=True)
    ax1.plot(time, 3*Fil_Ord, '--', color='black')
    ax1.plot(time, -3*Fil_Ord, '--', color='black')
    ax1.plot(time, res[:len(res),0], '-',color='red')
    ax1.fill_between(time, 3*Fil_Ord, -3* Fil_Ord, alpha=0.3)
    
    plt.title('Residual vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Residual')
    plt.show()
    
    
   





