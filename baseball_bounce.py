# Baseball Bounce Simulation
# For introducing a basic method of solving simple ODEs numerically
# REQUIRES PYTHON 3
# Uses blitting to speed up animation (doesn't work on Mac), see: https://matplotlib.org/3.3.0/tutorials/advanced/blitting.html

import numpy as np
from scipy import integrate
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import platform

# system parameters
sysParams = {
	"m": 0.145, # [kg]       total length of cable
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
 
	# speed up animation by skipping this many frames between refreshing plot
	anim_step = 5
 
	# integration time period
	t0 = 0.0
	tf = 1.2
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

	# prepare plot
	title_text = "Time: {:6.3f}s"
	ph_title = plt.title(title_text.format(0),fontweight='bold',animated=True)	
	plt.grid('on')
	plt.xlim(-0.309, 0.309)
	plt.ylim(-0.1,1)
	plt.gca().set_aspect('equal', adjustable='box')
	
	# draw background
	ground = patches.Rectangle((-0.309, -0.1),0.618,0.1,fill=True,facecolor='#999999',animated=False)
	plt.gca().add_patch(ground)

	# draw baseball
	(ph_ball,) = plt.gca().plot(0,np.nan,'o',markersize=20,markeredgewidth=5,markerfacecolor='#FFFFFF', markeredgecolor='#CC0000',animated=True)
	
	# draw and cache initial plot
	plt.show(block=False)
	plt.pause(0.1)
	bg = fig2.canvas.copy_from_bbox(fig2.bbox)
	plt.gca().draw_artist(ph_ball)
	plt.gca().draw_artist(ph_title)
	fig2.canvas.blit(fig2.bbox)

	# step through animation
	for i in np.arange(0,len(t.T),anim_step):
		
		# get current value of variable y (distance traveled)
		t_now = t.T[i][0]
		y = X[0,i]
	
		# update plot
		fig2.canvas.restore_region(bg)
		ph_ball.set_ydata(y)
		plt.gca().draw_artist(ph_ball)
		
		# add title
		ph_title.set_text(title_text.format(t_now))
		plt.gca().draw_artist(ph_title)
		
		# show current frame of animation
		fig2.canvas.blit(fig2.bbox)
		fig2.canvas.flush_events()
		#plt.pause(0.00001)          # don't plt.pause(), this is actually quite slow; instead change values of 'dt' and 'anim_step'
		
		# blitting doesn't work on Mac (known issue?) so if we're on a Mac add the pause
		if(platform.system() == "Darwin"):
			plt.pause(0.00001)		


	# keep figures displayed after animation ends
	plt.show(block=True)
