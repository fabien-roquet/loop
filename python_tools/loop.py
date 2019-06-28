
# This Python code implements the loop model of Roquet et al 2017

# It emulates the Matlab tools available at https://github.com/fabien-roquet/loop

 #Roquet, F., R. Lindqvist, F. Pollmann, D. Ferreira, and G. Madec (2017), 
 # Stability of the thermohaline circulation examined with a one-dimensional fluid loop, 
 #Tellus A: Dynamic Meteorology and Oceanography, 69(1), 
 #1380,490, doi: 10.1080/16000870.2017.1380490

# note that this script can read only output files in the format of fluiloop_v1 at the moment


import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve


def loop_equilibrium(w,R,nl,Zf,foldtrue,lamda,mu,ref_t,xi_t,F_s):
    # loop parameters    
    dl = 2*np.pi/nl
    l = ((np.array(range(nl))+1)*dl)
    
    # compute position of sink and source
    z = np.cos(l)
    n_sink = 0
    while (z[n_sink]-Zf)>= 1e-10 and n_sink<nl-1:
        n_sink = n_sink + 1
        
    n_source = nl - n_sink
    
    # compute mask
    mask = np.ones(nl)
    if foldtrue == 1:
        mask[0:n_sink] = 0
        mask[n_source-1: ] = 0
        
    # compute geometrical variables
    z = z*mask + (1-mask)*z[n_sink-1]
    curv = np.sin(l)*mask
    
    # compute analytically steady-state t/s distrib as a function of eccentric anomaly
    theta = loop_tracer_relax(w,R,nl,n_sink,ref_t,xi_t)
    salt = loop_tracer_fixed(w,R,nl,n_sink,F_s)
    sigma = -theta*(1+(lamda/2*theta).reshape(nl,1)-(mu*z).reshape(nl,1))+salt;
    torque = np.sum(sigma*curv.reshape(nl,1))/nl
    tau = w - torque
    return torque,tau,theta,salt

def loop_tracer_relax(w_0,R,nphi,n_sink,tracer_ref,xi):
    w=[]
    w.append(w_0)
    dphi = 2*np.pi/nphi
    phi = (np.array(range(nphi))+1)*dphi
    n_source = nphi - n_sink
    phi_sink = phi[n_sink-1]
    phi_source = phi[n_source-1]
    phi_diff = phi_source - phi_sink
    # avoid w = 0 case
    #if w[0] < 1e-10: w[0]=1e-10

    tracer = np.zeros([nphi,np.size(w)])
    
    for i in range(np.size(w)):
        # Solution for positive velocities
        if np.abs(w[i]) < 1e-10:
            w[i] = 1e-10
        if w[i] >= 1e-10:
            delt = R/w[i]
            phi_1 = np.where(np.logical_and(phi_sink<=phi,phi < phi_source))
            phi_2 = np.where(phi < phi_sink)
            phi_3 = np.where(phi>=phi_source)
            C = (1-np.exp((phi_diff-2*np.pi)/delt))/(1-np.exp(-phi_diff/delt))
            D = (1+np.exp((phi_diff-2*np.pi)/delt))/(1-np.exp((phi_diff-2*np.pi)/delt)+2*C*np.exp(-phi_diff/delt))
            t_inf_m = (2*np.pi*(xi/w[i])*tracer_ref/(-(1+D)*(np.exp((phi_diff-2*np.pi)/delt)+C)/(1+C*np.exp(-phi_diff/delt))
            -(2*np.pi*xi/w[i])*(D-np.exp((phi_diff-2*np.pi)/delt)*(1+D)/(1+C*np.exp(-phi_diff/delt)))))
            t_inf_p = -t_inf_m*D
            tracer[phi_1,i]=t_inf_m-(C*(t_inf_m-t_inf_p)/(1+C*np.exp(-phi_diff/delt)))*np.exp((phi[phi_1]-phi_source)/delt)
            tracer[phi_2,i]=t_inf_p+((t_inf_m-t_inf_p)/(1+C*np.exp(-phi_diff/delt)))*np.exp((phi[phi_2]-phi_sink)/delt)
            tracer[phi_3,i]=t_inf_p+((t_inf_m-t_inf_p)/(1+C*np.exp(-phi_diff/delt)))*np.exp((phi[phi_3]-phi_sink-2*np.pi)/delt)
        elif w[i] < 1e-10:
            delt = R/w[i]
            phi_1 = np.where(np.logical_and(phi_sink<=phi,phi < phi_source))
            phi_2 = np.where(phi < phi_sink)
            phi_3 = np.where(phi>=phi_source)
            C = (1-np.exp((2*np.pi-phi_diff)/delt))/(1-np.exp(phi_diff/delt))
            D = C*(1+np.exp(phi_diff/delt))/(1+np.exp((2*np.pi-phi_diff)/delt))
            t_inf_m = (2*np.pi*(xi/w[i])*tracer_ref/((1+D)*(1+C*np.exp(phi_diff/delt))/(1+C*np.exp(phi_diff/delt))
            -(2*np.pi*xi/w[i])*(D-C*np.exp(phi_diff/delt)*(1+D)/(1+C*np.exp(phi_diff/delt)))))
            t_inf_p = -t_inf_m*D
            tracer[phi_1,i]=t_inf_p-(C*(t_inf_p-t_inf_m)/(1+C*np.exp(phi_diff/delt)))*np.exp((phi[phi_1]-phi_sink)/delt)
            tracer[phi_2,i]=t_inf_m+((t_inf_p-t_inf_m)/(1+C*np.exp(phi_diff/delt)))*np.exp((phi[phi_2]+2*np.pi-phi_source)/delt)
            tracer[phi_3,i]=t_inf_m+((t_inf_p-t_inf_m)/(1+C*np.exp(phi_diff/delt)))*np.exp((phi[phi_3]-phi_source)/delt)
    return tracer

def loop_tracer_fixed(w_0,R,nphi,n_sink,flux):
    w=[]
    w.append(w_0)
    dphi = 2*np.pi/nphi
    phi = (np.array(range(nphi))+1)*dphi
    n_source = nphi - n_sink
    phi_sink = phi[n_sink-1]
    phi_source = phi[n_source-1]
    phi_diff = phi_source - phi_sink
    tracer = np.zeros([nphi,np.size(w)])
    
    for i in range(np.size(w)):
        # small velocity case
        if np.abs(w[i]) < 1e-5:
            w[i] = 1e-5
        # Solution for positive velocities
        if w[i] >= 1e-10:
            delt = R/w[i]
            phi_1 = np.where(np.logical_and(phi_sink<=phi,phi < phi_source))
            phi_2 = np.where(phi < phi_sink)
            phi_3 = np.where(phi>=phi_source)
            C = (1-np.exp((phi_diff-2*np.pi)/delt))/(1-np.exp(-phi_diff/delt))
            tracer[phi_1,i]=(flux*(phi_diff-2*np.pi)/w[i]+ 
                C*(flux*2*np.pi/w[i])*(np.exp((phi[phi_1]-phi_source)/delt)/(1+C*np.exp(-phi_diff/delt))))
            tracer[phi_2,i]=(flux*(phi_diff)/w[i]-
                (flux*2*np.pi/w[i])*(np.exp((phi[phi_2]-phi_sink)/delt)/(1+C*np.exp(-phi_diff/delt))))
            tracer[phi_3,i]=(flux*(phi_diff)/w[i]-
                (flux*2*np.pi/w[i])*(np.exp((phi[phi_3]-phi_sink-2*np.pi)/delt)/(1+C*np.exp(-phi_diff/delt))))
        # Solution for negative velocities
        elif w[i] < 1e-10:
            delt = R/w[i]
            phi_1 = np.where(np.logical_and(phi_sink<=phi,phi < phi_source))
            phi_2 = np.where(phi < phi_sink)
            phi_3 = np.where(phi>=phi_source)
            C = (1-np.exp(-(phi_diff-2*np.pi)/delt))/(1-np.exp(phi_diff/delt))
            tracer[phi_1,i]=(flux*(phi_diff-2*np.pi)/w[i]+
            C*(flux*2*np.pi/w[i])*(np.exp((phi[phi_1]-phi_sink)/delt)/(1+C*np.exp(phi_diff/delt))))
            tracer[phi_2,i]=(flux*(phi_diff)/w[i]-
            (flux*2*np.pi/w[i])*(np.exp((phi[phi_2]-phi_source+2*np.pi)/delt)/(1+C*np.exp(phi_diff/delt))))
            tracer[phi_3,i]=(flux*(phi_diff)/w[i]-
            (flux*2*np.pi/w[i])*(np.exp((phi[phi_3]-phi_source)/delt)/(1+C*np.exp(phi_diff/delt))))

    return tracer
    
def loop_jacobian(w,R,nl,Zf,foldtrue,lamda,mu,ref_t,xi_t,F_s):
    dl = 2*np.pi/nl
    l = (np.array(range(nl))+1)*dl
    z = np.cos(l)
    n_sink = 0
    while (z[n_sink]-Zf) >= 1e-10 and n_sink < nl-1:
        n_sink = n_sink+1
    n_source = nl - n_sink
    
    # compute mask
    mask = np.ones(nl)
    if foldtrue == 1:
        mask[0:n_sink] = 0
        mask[n_source-1: ] = 0
        
    # compute geometrical variables
    z = z*mask + (1-mask)*z[n_sink-1]
    curv = (np.sin(l)*mask).reshape(nl,1)
        
    # compute analytically steady-state t/s distrib as a function of eccentric anomaly
    theta = loop_tracer_relax(w,R,nl,n_sink,ref_t,xi_t)
    salt = loop_tracer_fixed(w,R,nl,n_sink,F_s)
    
    # Compute Jacobian
    # init matrix 
    J = np.zeros([2*nl,2*nl])
    
    # temperature
    for nn in range(nl): J[nn,nn] = -2*R/dl**2
    for nn in range(nl-1): J[nn,nn+1]= (w/2 + R/dl)/dl
    J[nl-1,0] = (w/2 + R/dl)/dl
    for nn in range(nl-1): J[nn+1,nn] = (-w/2 + R/dl)/dl
    J[0,nl-1] =  (-w/2 + R/dl)/dl
    J[n_source-1,n_source-1] = J[n_source-1,n_source-1] - ref_t*xi_t*nl
    J[n_sink-1,n_sink-1] = J[n_sink-1,n_sink-1] - ref_t*xi_t*nl
    for nm in range(nl):
        J[nm,0] = J[nm,0]+curv[nm]*(1+lamda*theta[nm]-mu*z[nm])/nl*(theta[1]-theta[nl-1])/(2*dl)
        J[nm+nl,0] = J[nm+nl,0]-curv[nm]/nl*(theta[1]-theta[nl-1])/(2*dl)
        for nn in range(1,nl-1):
            J[nm,nn]=J[nm,nn]+curv[nm]*(1+lamda*theta[nm]-mu*z[nm])/nl*(theta[nn+1]-theta[nn-1])/(2*dl)
            J[nm+nl,nn]=J[nm+nl,nn]-curv[nm]/nl*(theta[nn+1]-theta[nn-1])/(2*dl)
        J[nm,nl-1]=J[nm,nl-1]+curv[nm]*(1+lamda*theta[nm]-mu*z[nm])/nl*(theta[0]-theta[nl-2])/(2*dl)
        J[nm+nl,nl-1]=J[nm+nl,nl-1]-curv[nm]/nl*(theta[0]-theta[nl-2])/(2*dl)
        
    # salinity
    for nn in range(nl,2*nl): J[nn,nn]= -2*R/dl**2
    for nn in range(nl+1,2*nl): J[nn-1,nn]= (w/2+R/dl)/dl
    J[2*nl-1,nl]=(w/2+R/dl)/dl
    for nn in range(nl,2*nl-1): J[nn+1,nn]= (-w/2+R/dl)/dl
    J[nl,2*nl-1]=(-w/2+R/dl)/dl
    for nm in range(nl):
        J[nm,nl]=J[nm,nl]+curv[nm]*(1+lamda*theta[nm]-mu*z[nm])/nl*(salt[1]-salt[nl-1])/(2*dl)
        J[nm+nl,nl]=J[nm+nl,nl]-curv[nm]/nl*(salt[1]-salt[nl-1])/(2*dl)
        for nn in range(1,nl-1):
            J[nm,nn+nl]=J[nm,nn+nl]+curv[nm]*(1+lamda*theta[nm]-mu*z[nm])/nl*(salt[nn+1]-salt[nn-1])/(2*dl)
            J[nm+nl,nn+nl]=J[nm+nl,nn+nl]-curv[nm]/nl*(salt[nn+1]-salt[nn-1])/(2*dl)
        J[nm,2*nl-1]=J[nm,2*nl-1]+curv[nm]*(1+lamda*theta[nm]-mu*z[nm])/nl*(salt[0]-salt[nl-2])/(2*dl)
        J[nm+nl,2*nl-1]=J[nm+nl,2*nl-1]-curv[nm]/nl*(salt[0]-salt[nl-2])/(2*dl)
    return J
    
def loop_init_state(file_init_state,w,theta,salt):
    nl = len(theta)
    fid = open(file_init_state,'ab')    
    np.savetxt(fid,np.column_stack((np.arange(1,nl+1).astype(int),theta.reshape(nl,),salt.reshape(nl,))),fmt='%i %2.5f %2.5f')
    fid.close()
    
def loop_read_out(name_exp):
    name_out_0d = '%s_out_0D.txt'%name_exp
    name_out_1d = '%s_out_1D.txt'%name_exp
    out_0d = {}
    out_1d = {}    
    
    with open(name_out_0d,'r') as f0:        
        y0=[[float(s) for s in line.split( )] for line in f0.readlines()[1:]]
    out_0d['niter'] = np.array(y0)[:,0]
    out_0d['time'] = np.array(y0)[:,1]
    out_0d['w'] = np.array(y0)[:,2]
    out_0d['mass'] = np.array(y0)[:,3]

    with open(name_out_1d,'r') as f1:        
        y1=[[float(s) for s in line.split( )] for line in f1.readlines()[1:]]
    niter = np.array(y1)[:,0]
    jk = np.array(y1)[:,1]
    theta = np.array(y1)[:,2]
    salt = np.array(y1)[:,3]
    sigma = np.array(y1)[:,4]
    nl = int(max(jk))
    nk = int(len(niter)/nl)
    out_1d['niter'] = niter.reshape(nk,nl)
    out_1d['theta'] = theta.reshape(nk,nl)
    out_1d['salt'] = salt.reshape(nk,nl)
    out_1d['sigma'] = sigma.reshape(nk,nl)
    return out_0d,out_1d
    
# test configuration    
nl = 360
Zf = .5
foldtrue = 1
R = .1
lamda = 0
mu = 0
ref_t = 5
xi_t = 1
F_s = .15
tau = 0

# plot torque curves
w_list = np.arange(-1,1.01,0.01)
torque = w_list*np.nan
for kk in range(np.size(w_list)):
    torque[kk] = loop_equilibrium(w_list[kk],R,nl,Zf,foldtrue,lamda,mu,ref_t,xi_t,F_s)[0]
    
# determine semi-analytically steady-state
func = lambda w : (w - loop_equilibrium(w,R,nl,Zf,foldtrue,lamda,mu,ref_t,xi_t,F_s)[0] - tau)
w1 = fsolve(func,-.1)
torque1 = w1 - tau;
w2 = fsolve(func, .3)
torque2 = w2 - tau;
w3 = fsolve(func, .1)
torque3 = w3 - tau
  
# determine the jacobian and its eigenvalues for steady-states 1
jacobian1 = loop_jacobian(w1,R,nl,Zf,foldtrue,lamda,mu,ref_t,xi_t,F_s)  
EIG1=np.linalg.eig(jacobian1)[0]
lEIG1=EIG1[np.where(EIG1.real > 0)]
print('STATE 1, w =%4.2f : no eigenvalues with positive real part --> stable state'%w1)

# determine the jacobian and its eigenvalues for steady-states 2
jacobian2 = loop_jacobian(w2,R,nl,Zf,foldtrue,lamda,mu,ref_t,xi_t,F_s)  
EIG2=np.linalg.eig(jacobian2)[0]
lEIG2=EIG2[np.where(EIG2.real > 0)]
print('STATE 2, w =%4.2f : two conjugate eigenvalues with positive real part %4.2f +/- %4.2fi --> oscillatory state'%(w2,lEIG2[0].real,lEIG2[0].imag))

# determine the jacobian and its eigenvalues for steady-states 3
jacobian3 = loop_jacobian(w3,R,nl,Zf,foldtrue,lamda,mu,ref_t,xi_t,F_s)  
EIG3=np.linalg.eig(jacobian3)[0]
lEIG3=EIG3[np.where(EIG3.real > 0)]
print('STATE 3, w =%4.2f : one real positive eigenvalue %4.2f --> unstable state'%(w3,lEIG3.real))

#plot torque curves. Each intersection is a steady-state
fig1,c1 = plt.subplots()
c1.plot(w_list[:-1],w_list[:-1],label=r'$\tau$')
c1.plot(w_list[:-1],torque[:-1],label='buoyancy torque')
c1.plot(w1,w1,'ro')
c1.plot(w2,w2,'ro')
c1.plot(w3,w3,'ro')
c1.set_xlabel('velocity')
c1.set_ylabel('torque')
c1.legend(loc='lower right')
plt.title('Torque curves')
plt.show()

[torque1,tau1,theta1,salt1] = loop_equilibrium(w1,R,nl,Zf,foldtrue,lamda,mu,ref_t,xi_t,F_s)
loop_init_state('./init_state.txt',w1,theta1,salt1)

[out_0d,out_1d]=loop_read_out('TEST')
print('\n')
print('Difference between loop output and theoretical solution of STATE1 velocity   : %4.2g'%(out_0d['w'][-1]-w1))
print('STD between loop output and theoretical solution of STATE1 theta distribution: %4.2g'%np.std(out_1d['theta'][-1,:]-np.transpose(theta1)))
print('STD between loop output and theoretical solution of STATE1 salt distribution: %4.2g'%np.std(out_1d['salt'][-1,:]-np.transpose(salt1)))

fig2,c2 = plt.subplots()
xaxis = np.arange(1,nl+1)
c2.plot(xaxis,out_1d['theta'][-1,:],label='STATE1 theta')
c2.plot(out_1d['salt'][-1,:],label='STATE1 salt')
c2.legend(loc='lower right')
plt.show()
