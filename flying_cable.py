# Flying Cable Simulation
# Feynman exercise #10.15
# REQUIRES PYTHON 3
# Uses blitting to speed up animation (doesn't work on Mac), see: https://matplotlib.org/3.3.0/tutorials/advanced/blitting.html

import numpy as np
from scipy import integrate
import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import platform

# system parameters
sysParams = {
	"L": 1, 	# [m]       total length of cable
	"g": 9.81,  # [m/s^2]   acceleration of gravity
	"rho": 1,   # [kg/m]    rope mass per unit length
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
	r = integrate.ode(StateProp).set_integrator('dopri5')
 
 	# speed up animation by skipping this many frames between refreshing plot
	anim_step = 5
 
	# integration time period
	t0 = 0.0
	tf = 1.5
	dt = 0.001
	nSteps = int(np.floor(((tf-t0)/dt)+1))
	
	# initial conditions
	X0 = np.zeros((2,1))
	X0[0] = 0.005 # [m]
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
	
	# Mac OS is silly and doesn't support blitting
	# so run matplotlib with the TkAgg backend
	# see: https://retifrav.github.io/blog/2020/09/05/matplotlib-animation-macos/
	if(platform.system() == "Darwin"):
		matplotlib.use("TkAgg")
	
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
	fig.tight_layout()
	fig.canvas.draw() 
		
	# animate results
	fig2 = plt.figure(num=None,figsize=(4,10))
	r_pulley = 0.1  # for display purposes only... in developing the physics we assumed a negligible pulley radius 
	
	# prepare plot
	title_text = "Time: {:6.3f}s"
	ph_title = plt.title(title_text.format(0),fontweight='bold',animated=True)
	plt.grid('on')
	plt.xlim(-0.4, 0.4)
	plt.ylim(-2.5,2*r_pulley)
	plt.gca().set_aspect('equal', adjustable='box')

	# draw background
	theta = np.arange(0,np.pi,0.01)
	x_circ = r_pulley*np.cos(theta)
	y_circ = r_pulley*np.sin(theta)
	pulley = patches.Circle((0, 0), 0.8*r_pulley,fill=True,facecolor='#000000',animated=False)
	plt.gca().add_patch(pulley)
	
	# draw three segments of cable
	seg_nan_lr = np.empty([2,1])
	seg_nan_lr[:] = np.nan
	seg_nan_ctr = np.empty(y_circ.shape)
	seg_nan_ctr[:] = np.nan
	(ph_lft_seg,) = plt.plot(-1*r_pulley*np.ones([2,1]),seg_nan_lr,'-',linewidth=5,color='#cc0000',animated=True)
	(ph_ctr_seg,) = plt.plot(x_circ,seg_nan_ctr,'-',linewidth=5,color='#cc0000',animated=True)
	(ph_rht_seg,) = plt.plot(r_pulley*np.ones([2,1]),seg_nan_lr,'-',linewidth=5,color='#cc0000',animated=True)
	
	# draw and cache initial plot
	plt.show(block=False)
	plt.pause(0.1)
	bg = fig2.canvas.copy_from_bbox(fig2.bbox)
	plt.gca().draw_artist(ph_lft_seg)
	plt.gca().draw_artist(ph_ctr_seg)
	plt.gca().draw_artist(ph_rht_seg)
	plt.gca().draw_artist(ph_title)
	fig2.canvas.blit(fig2.bbox)
	
	# step through animation
	for i in np.arange(0,len(t.T),anim_step):
		
		# get current value of variable y (distance traveled)
		t_now = t.T[i][0]
		y = X[0,i]
	
		# update plot
		fig2.canvas.restore_region(bg)
		
		# show cable
		if(y <= sysParams["L"]/2):
			ph_lft_seg.set_ydata(np.array([[y-(sysParams["L"]/2)],[0]]))
			ph_ctr_seg.set_ydata(y_circ)
			ph_rht_seg.set_ydata(np.array([[0],[-(sysParams["L"]/2+y)]]))
			plt.gca().draw_artist(ph_lft_seg)
			plt.gca().draw_artist(ph_ctr_seg)
			plt.gca().draw_artist(ph_rht_seg)			
		else:
			ph_rht_seg.set_ydata((-y+(sysParams["L"]/2))+np.array([[0],[-sysParams["L"]]]))
			plt.gca().draw_artist(ph_rht_seg)
			
		# add title
		ph_title.set_text(title_text.format(t_now))
		plt.gca().draw_artist(ph_title)
		
		# show current frame of animation
		fig2.canvas.blit(fig2.bbox)
		fig2.canvas.flush_events()
		#plt.pause(0.00001)          # actually don't use plt.pause(), this is actually quite slow; instead change values of 'dt' and 'anim_step'

	# keep figures displayed after the animation ends
	# unfortunately plt.show(block=True) clears the blitted animation in Mac OS with TkAgg backend
	# so we'll just wait for a keypress instead
	print("Press any key to exit.")
	input()
