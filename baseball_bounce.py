# Baseball Bounce Simulation
# For introducing basics of using ode45 for solving simple ODEs numerically
# REQUIRES PYTHON 3

import numpy as np
from scipy import integrate
import matplotlib.patches as patches
import matplotlib.pyplot as plt

# system parameters
sysParams = {
	"m": 1, 	# [kg]       total length of cable
	"g": 9.81,  # [m/s^2]    acceleration of gravity
	"e": 0.45,  # coefficient of restitution
}
 
# derivative calculation
# this function is called by the ODE solver
def StateProp(t, X):
	"""
	Calculate Xdot at any time t
	"""
	
	# deconstruct state vector
	y = X[0]
	y_dot = X[1]
	
	# construct Xdot from differential equation
	# note: X = [y y_dot] therefore Xdot = [y_dot y_ddot]
	Xdot = np.zeros((2,1))
	Xdot[0] = y_dot	
	Xdot[1] = -1*sysParams["g"]
	
	return Xdot

 
# main function for running ODE solver
if __name__ == '__main__':
 
	# integrator selection
	# dopri5 is a Dormand-Price Runge-Kutta 4/5 solver
	r = integrate.ode(StateProp).set_integrator('dopri5')
 
	# integration time period
	t0 = 0.0
	tf = 1.15
	dt = 0.001
	nSteps = int(np.floor(((tf-t0)/dt)+1))
	
	# initial conditions
	X0 = np.zeros((2,1))
	X0[0] = 1.0   # [m]
	X0[1] = 0     # [m/s]
	r.set_initial_value(X0, t0)
 
	# data storage
	t = np.zeros((1,nSteps))
	X = np.zeros((2,nSteps))
	t[0] = 0
	X[:,0] = X0.T    # silly numpy matrix style requires that column to be inserted is actually a row...
 
	# integrate for each timestep and store output
	k = 1
	while r.successful() and k < nSteps:
	
		# apply bounce if needed
		if( r.y[1] < 0 and r.y[0] < 0):
			r.y[1] = -1*sysParams["e"]*r.y[1]
	
		r.integrate(r.t + dt)
		t[0,k] = r.t
		X[:,k] = r.y.T
		k += 1
 

	# compute total energy in system (i.e. the Hamiltonian)
	# depends upon whether the cable is unwrapping around pulley (E1)
	# or falling freely (E2)
	y = X[0,:]
	y_dot = X[1,:]
	E = 0.5*sysParams["m"]*y_dot**2 + sysParams["m"]*sysParams["g"]*y;
	
	# plot results using matplotlib
	t_bounds = [[t.T[0][0]],[t.T[len(t.T)-1][0]]]
	fig = plt.figure(num=None,figsize=(10, 8))
	plt.subplot(311)
	plt.plot(t.T, X[0,:],'-',linewidth=3,color='#0000CC')
	plt.grid('on')
	plt.xlabel('Time [s]',fontweight='bold')
	plt.ylabel('Position [m]',fontweight='bold')
	plt.subplot(312)
	plt.plot(t.T, X[1,:],'-',linewidth=3,color='#CC0000')
	plt.grid('on')
	plt.xlabel('Time [s]',fontweight='bold')
	plt.ylabel('Speed [m/s]',fontweight='bold')	
	plt.show(block=False)
	plt.subplot(313)
	plt.plot(t.T,E,'-',linewidth=3,color='#00CC00')
	plt.grid('on')
	plt.xlabel('Time [s]',fontweight='bold')
	plt.ylabel('Energy [J]',fontweight='bold')	
	fig.tight_layout()
	fig.canvas.draw() 
		
	# animate results
	fig2 = plt.figure(num=None,figsize=(4,10))

	for i in np.arange(0,len(t.T),10):
		
		# get current value of variable y (distance traveled)
		t_now = t.T[i][0]
		y = X[0,i]
	
		# format plot
		fig2.clear()
		title_text = "Time: {:6.3f}s"
		plt.title(title_text.format(t_now),fontweight='bold')
		plt.grid('on')
		plt.xlim(-0.309, 0.309)
		plt.ylim(-0.1,1)
		plt.gca().set_aspect('equal', adjustable='box')
		
		# show ground
		art = patches.Rectangle((-0.309, -0.1),0.618,0.1,fill=True,facecolor='#999999')
		plt.gca().add_patch(art)

		# show baseball`
		plt.plot(0,y,'o',markersize=20,markeredgewidth=5,markerfacecolor='#FFFFFF', markeredgecolor='#CC0000')
		
		# show current frame of animation
		fig2.canvas.draw()
		fig2.canvas.flush_events()
		plt.pause(0.0001)
	
	# keep figures displayed after animation ends
	plt.show(block=True)
