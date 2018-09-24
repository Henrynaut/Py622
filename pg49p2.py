#Written by Neil McHenry
#September 23rd, 2018

#Load math relate packages
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import ode

#AERO 622 Hw #3 Page 49 Problem 2
    #Solve for the motion of the bead on a spinning wire
    #Using a numerical integration technique, then comment

#Physical parameters
m = 0.1     #kg
r = 0.4     #meters
g = 10      #m/s/s

#Initial Condition 1
OMEGA1 = 4               #rad/s
theta = np.pi/2 + 0.1 #radians
thetaDot = 0            #rad/s
time = 12               #seconds

#Initial Condition 2
OMEGA2 = 6   #rad/s

#Function Definition
def fun(t, z, omega):
    """
    Right hand side of the differential equations
      dx/dt = -omega * y
      dy/dt = omega * x
    """
    x, y = z
    f = [-omega*y, omega*x]
    return f

# Create an `ode` instance to solve the system of differential
# equations defined by `fun`, and set the solver method to 'dop853'.
solver = ode(fun)
solver.set_integrator('dop853')

# Give the value of omega to the solver. This is passed to
# `fun` when the solver calls it.
omega = 2 * np.pi
solver.set_f_params(omega)

# Set the initial value z(0) = z0.
t0 = 0.0
z0 = [1, -0.25]
solver.set_initial_value(z0, t0)

# Create the array `t` of time values at which to compute
# the solution, and create an array to hold the solution.
# Put the initial value in the solution array.
t1 = 2.5
N = 75
t = np.linspace(t0, t1, N)
sol = np.empty((N, 2))
sol[0] = z0

# Repeatedly call the `integrate` method to advance the
# solution to time t[k], and save the solution in sol[k].
k = 1
while solver.successful() and solver.t < t1:
    solver.integrate(t[k])
    sol[k] = solver.y
    k += 1

# Plot the solution...
plt.plot(t, sol[:,0], label='x')
plt.plot(t, sol[:,1], label='y')
plt.xlabel('t')
plt.grid(True)
plt.legend()
plt.show()