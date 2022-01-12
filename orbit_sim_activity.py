# Simple Orbit Simulation
# Author: Mike Kokko
# Modified: 11-Jan-2022
#
# For introducing a basic method of solving simple ODEs numerically
# REQUIRES PYTHON 3
# Uses blitting to speed up animation (probably doesn't work on Mac), see: https://matplotlib.org/3.3.0/tutorials/advanced/blitting.html
#
import numpy as np
from scipy import integrate
import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import platform
import math
import numpy

# system parameters
sysParams = {
    "M": 5.976e24, # [kg]            mass of Earth
    "G": 6.67e-11, # [N*m^2/kg^2]    gravitational constant
    "re": 6.373e6, # [m]             radius of Earth
    "done": False
}
 
# derivative calculation
# this function is called by the ODE solver
def StateProp(t, X):
    """
    Calculate Xdot at any time t
    """
    
    # deconstruct state vector
    r = X[0]
    r_dot = X[1]
    theta_dot = X[3]
    
    # construct Xdot from differential equation
    # note:     X    = [y y_dot]
    # therefore Xdot = [y_dot y_ddot]
    Xdot = np.zeros((4,1))
    Xdot[0] = r_dot
    Xdot[1] = <<ENTER R DOUBLE DOT EXPRESSION HERE>>
    Xdot[2] = theta_dot
    Xdot[3] = <<ENTER THETA DOUBLE DOT EXPRESSION HERE>> 
    
    return Xdot

 
# main function for running ODE solver
if __name__ == '__main__':

    # integrator selection
    # dopri5 is a Dormand-Price Runge-Kutta 4/5 solver
    r = integrate.ode(StateProp).set_integrator('dopri5')
 
    # speed up animation by skipping this many frames between refreshing plot
    anim_step = 20
 
    # integration time period
    t0 = 0.0
    tf = 200000.0
    dt = 5
    nSteps = int(np.floor(((tf-t0)/dt)+1))
    
    # initial conditions
    r_0 = <<ENTER RADIUS HERE>> # [m]
    r_dot_0 = 0 # [m/s]    
    theta_0 = 0 # [rad]
    theta_dot_0 = <<ENTER SPEED HERE>>/(r_0) # [rad/s]

    # apply initial conditions
    X0 = np.zeros((4,1))
    X0[0] = r_0          # [m]
    X0[1] = r_dot_0      # [m/s]
    X0[2] = theta_0      # [rad]
    X0[3] = theta_dot_0  # [rad/s]    
    r.set_initial_value(X0, t0)
 
    # data storage
    t = np.zeros((1,nSteps))
    X = np.zeros((4,nSteps))
    t[0] = 0
    X[:,0] = X0.T    # silly numpy matrix style requires that column to be inserted is actually a row...
 
    # integrate for each timestep and store output
    k = 1
    while r.successful() and k < nSteps and (not sysParams["done"]):
    
        # check for collision with Earth or complete orbit
        if((X[0,k-1] < sysParams["re"]) or (X[2,k-1] > 2*math.pi)):
            sysParams["done"] = True
        
        # propagate state
        r.integrate(r.t + dt)
        t[0,k] = r.t
        X[:,k] = r.y.T
        k += 1
 

    # TODO: compute total energy in system (i.e. the Hamiltonian)
    
    # Mac OS is silly and doesn't support blitting
    # so run matplotlib with the TkAgg backend
    # see: https://retifrav.github.io/blog/2020/09/05/matplotlib-animation-macos/
    if(platform.system() == "Darwin"):
        matplotlib.use("TkAgg")

    # extract and compute orbit trajectory parameters
    r = X[0,0:k-1]
    r_dot = X[1,0:k-1]
    theta = X[2,0:k-1]
    theta_dot = X[3,0:k-1]
    x = r*numpy.cos(theta)
    y = r*numpy.sin(theta)
    v = numpy.sqrt( (r_dot)**2 + (r*theta_dot)**2 )
 
    # plot position using matplotlib
    t_bounds = [[t.T[0][0]],[t.T[len(t.T)-1][0]]]
    fig = plt.figure(num=None)
    plt.subplot(211)
    plt.plot(t.T[0:k-1]/60, r,'-',linewidth=3,color='#0000CC')
    plt.grid('on')
    plt.xlabel('Time [min]',fontweight='bold')
    plt.ylabel('Radius [m]',fontweight='bold')

    # plot speed using matplotlib    
    plt.subplot(212)
    plt.plot(t.T[0:k-1]/60, v,'-',linewidth=3,color='#0000CC')
    plt.grid('on')
    plt.xlabel('Time [min]',fontweight='bold')
    plt.ylabel('Speed [m/s]',fontweight='bold')
    
    # finish formatting plot
    fig.tight_layout()
    fig.canvas.draw() 
        
    # animate results
    fig2 = plt.figure(num=None,figsize=(10,6))
    
    # prepare plot
    title_text = "Time: {:6.0f}min"
    ph_title = plt.title(title_text.format(0),fontweight='bold',animated=True)    
    plt.grid('on')
    #plt.xlim(-3e7, 3e7)
    #plt.ylim(-e3,1)
    ax = plt.gca()
    ax.set(facecolor = "black")
    plt.gca().set_aspect('equal', adjustable='box')
    
    # draw background
    earth = patches.Circle((0,0),radius=sysParams["re"],fill=True,facecolor='#6666FF',animated=False)
    plt.gca().add_patch(earth)
    plt.plot(x, y,'-',linewidth=3,color='#FFFF00')
 
    # draw satellite
    plt.gca().plot(x[0],y[0],'o',markersize=10,markeredgewidth=5,markerfacecolor='#FF00FF', markeredgecolor='#FF00FF',animated=False)
    (ph_sat,) = plt.gca().plot(0,np.nan,'o',markersize=5,markeredgewidth=5,markerfacecolor='#FFFF00', markeredgecolor='#FFFF00',animated=True)
    
    # draw and cache initial plot
    plt.show(block=False)
    plt.pause(0.1)
    bg = fig2.canvas.copy_from_bbox(fig2.bbox)
    plt.gca().draw_artist(ph_sat)
    plt.gca().draw_artist(ph_title)
    fig2.canvas.blit(fig2.bbox)

    # step through animation
    for i in np.arange(0,k-1,anim_step):
        
        # get current value of variable y (distance traveled)
        t_now = t.T[i][0]/60
        sat_x = x[i]
        sat_y = y[i]
    
        # update plot
        fig2.canvas.restore_region(bg)
        ph_sat.set_xdata(sat_x)
        ph_sat.set_ydata(sat_y)
        plt.gca().draw_artist(ph_sat)
        
        # add title
        ph_title.set_text(title_text.format(t_now))
        plt.gca().draw_artist(ph_title)
        
        # show current frame of animation
        fig2.canvas.blit(fig2.bbox)
        fig2.canvas.flush_events()
        #plt.pause(0.00001)          # on second thought, don't use plt.pause(), this is actually quite slow; instead change values of 'dt' and 'anim_step'
        
    # keep figures displayed after the animation ends
    # unfortunately plt.show(block=True) clears the blitted animation in Mac OS with TkAgg backend
    # so we'll just wait for a keypress instead
    print("Press any key to exit.")
    input()