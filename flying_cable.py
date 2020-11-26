# drawn largely from http://modelling3e4.connectmv.com/wiki/Software_tutorial/Integration_of_ODEs

import numpy as np
from scipy import integrate
import matplotlib.patches as patches
import matplotlib.pyplot as plt

# system parameters
sysParams = {
	"L": 1, 	# [m]       total length of cable
	"g": 9.81,  # [m/s^2]   acceleration of gravity
	"rho": 1,   # [kg/m]    rope mass per unit length
}
 
# derivative calculation
# this function is called by the ODE solver
def pendulumTestStateProp(t, X):
	"""
	Calculate Xdot at any time t
	"""
	
	# deconstruct state vector
	y = X[0]
	y_dot = X[1]
	
	# construct Xdot from differential equation
	# note: X = [y y_dot] therefore Xdot = [y_dot y_ddot]
	Xdot = np.zeros((2,1))
	Xdot[0] = y_dot;
	
	if(y <= sysParams["L"]/2):
		Xdot[1] = 2*y*sysParams["g"]/sysParams["L"]
	else:
		Xdot[1] =  sysParams["g"]

	return Xdot
 
# main function for running ODE solver
if __name__ == '__main__':
 
	# integrator selection
	# dopri5 is a Dormand-Price Runge-Kutta 4/5 solver
	r = integrate.ode(pendulumTestStateProp).set_integrator('dopri5')
 
	# integration time period
	t0 = 0.0
	tf = 1.5
	dt = 0.001
	nSteps = int(np.floor(((tf-t0)/dt)+1))
	
	# initial conditions
	X0 = np.zeros((2,1))
	X0[0] = 0.005 # degrees
	X0[1] = 0     # degrees/sec
	r.set_initial_value(X0, t0)
 
	# data storage
	t = np.zeros((1,nSteps))
	X = np.zeros((2,nSteps))
	t[0] = 0
	X[:,0] = X0.T    # silly numpy matrix style requires that column to be inserted is actually a row...
 
	# integrate for each timestep and store output
	k = 1
	while r.successful() and k < nSteps:
		r.integrate(r.t + dt)
		t[0,k] = r.t
		X[:,k] = r.y.T
		k += 1
 
	# figure out when cable leaves pulley
	offPulleyIdx = np.where(X[0,:] >= (sysParams["L"]/2))
	targetIdx = offPulleyIdx[0][0]	
	txt = "Cable leaves pulley at t = {:5.3f}s with v = {:5.3f}m/s"
	print(txt.format(t.T[targetIdx][0],X[1,targetIdx]))

	# compute total energy in system (i.e. the Hamiltonian)
	# depends upon whether the cable is unwrapping around pulley (E1)
	# or falling freely (E2)
	y = X[0,:]
	y_dot = X[1,:]
	m1 = (sysParams["L"]/2)-y
	m2 = (sysParams["L"]/2)+y
	m = m1+m2 # (== sysParams["L"] !)
	E1 = 0.5*(sysParams["rho"]*sysParams["L"])*np.square(y_dot) + sysParams["rho"]*sysParams["g"]*(0.25*sysParams["L"]**2-np.square(y));
	E2 = 0.5*(sysParams["rho"]*sysParams["L"])*np.square(y_dot) + sysParams["rho"]*sysParams["g"]*(0.5*sysParams["L"]-y);
	E1_idx = offPulleyIdx = np.where(X[0,:] < (sysParams["L"]/2))
	E2_idx = offPulleyIdx = np.where(X[0,:] >= (sysParams["L"]/2))
	E = np.concatenate((E1[E1_idx],E2[E2_idx]))
	
	# plot results using matplotlib
	t_bounds = [[t.T[0][0]],[t.T[len(t.T)-1][0]]]
	fig = plt.figure(num=None,figsize=(10, 8))
	plt.subplot(311)
	plt.plot(t_bounds,X[0,targetIdx]*np.ones([2,1]),'-',linewidth=3,color='#CC00CC')
	plt.plot(t.T, X[0,:],'-',linewidth=3,color='#0000CC')
	plt.grid('on')
	plt.xlabel('Time [s]',fontweight='bold')
	plt.ylabel('Position [m]',fontweight='bold')
	plt.subplot(312)
	plt.plot(t_bounds,X[1,targetIdx]*np.ones([2,1]),'-',linewidth=3,color='#CC00CC')
	plt.plot(t.T, X[1,:],'-',linewidth=3,color='#CC0000')
	plt.grid('on')
	plt.xlabel('Time [s]',fontweight='bold')
	plt.ylabel('Velocity [m/s]',fontweight='bold')	
	plt.show(block=False)
	plt.subplot(313)
	plt.plot(t.T,E,'-',linewidth=3,color='#00CC00')
	plt.grid('on')
	plt.xlabel('Time [s]',fontweight='bold')
	plt.ylabel('Energy [J]',fontweight='bold')	
	
	fig.canvas.draw() 
	

	
	# animate results
	r_pulley = 0.1
	theta = np.arange(0,np.pi,0.01)
	x_circ = r_pulley*np.cos(theta)
	y_circ = r_pulley*np.sin(theta)
	fig2 = plt.figure(num=None,figsize=(4,10))

	for i in np.arange(0,len(t.T),20):
		
		# get current value of variable y (distance traveled)
		t_now = t.T[i][0]
		y = X[0,i]
	
		# format plot
		fig2.clear()
		title_text = "Time: {:6.3f}s"
		plt.title(title_text.format(t_now),fontweight='bold')
		plt.grid('on')
		plt.xlim(-0.4, 0.4)
		plt.ylim(-2.5,2*r_pulley)
		plt.gca().set_aspect('equal', adjustable='box')
		
		# show pulley
		art = patches.Circle((0, 0), 0.8*r_pulley,fill=True,facecolor='#000000')
		plt.gca().add_patch(art)
		
		# show cable
		if(y <= sysParams["L"]/2):
			plt.plot(-1*r_pulley*np.ones([2,1]),np.array([[y-(sysParams["L"]/2)],[0]]),'-',linewidth=5,color='#cc0000');
			plt.plot(x_circ,y_circ,'-',linewidth=5,color='#cc0000');
			plt.plot(r_pulley*np.ones([2,1]),np.array([[0],[-(sysParams["L"]/2+y)]]),'-',linewidth=5,color='#cc0000')
		else:
			plt.plot(r_pulley*np.ones([2,1]),(-y+(sysParams["L"]/2))+np.array([[0],[-sysParams["L"]]]),'-',linewidth=5,color='#cc0000')
		
		# show current frame of animation
		fig2.canvas.draw()
		fig2.canvas.flush_events()
		plt.pause(0.0001)
	
	# keep figures displayed after animation ends
	plt.show(block=True)
