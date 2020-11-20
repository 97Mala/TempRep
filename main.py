import numpy as np
from scipy.integrate import ode
import os
#  set the 'MPLCONFIGDIR' to '/tmp' to get a writable environment for matplotlib
os.environ['MPLCONFIGDIR'] = '/tmp'
import matplotlib.pyplot as plt
from KINmodule import *

# time discretizzation
t0 = 0
tf = 60
dt = 0.1 # [s]
N   = round((tf-t0)/dt)+1
t = np.array(range(N))*dt+t0

# data loading
[u0, L, l, b] = dataBuilder()

# initializing matrix A: statedot system matrix (input dependent)
tmp = [[0.0 for i in range(7)] for j in range(7)]  # temporary array
A = np.array(tmp)
# constant part of matrix A
for i in range(6):
    A[(i + 1), 0] = b[i] / L
    A[(i + 1), (i + 1)] = -l[i]
# set the STATEDOT parameters (see KINmodule)
# REACTIVITY FUNCTION CHOICE HERE (first argument of param)
param = [rho_step, L, l, b,A]

######################################################################
# problem solution ###################################################

# SOLUTION initializzation
ts = []
sol = []
# SOLUTOR parameters (same solutor as MATLAB ode15s )
ode15s = ode(statedot)
# max step limitation necessary
ode15s.set_integrator('zvode', method='bdf', order=15, max_step=dt / 7)
# setting initial values and reactivity input function
ode15s.set_initial_value(u0, t=0.0).set_f_params(param)
# (DN is the number of steps between two message)
DN = 100
nsteps = DN
print('------------------------------')
# SOLVE THE PROBLEM
while ode15s.successful() and ode15s.t < tf:
    # this tells us how many steps have been performed
    if (i / nsteps == 1):
        print('------------------------------')
        print(str(nsteps), ' steps performed')
        nsteps = nsteps + DN
        print('------------------------------')
    # stepwise
    ode15s.integrate(ode15s.t + dt)
    # loading
    ts.append(ode15s.t)
    sol.append(np.real(ode15s.y.copy()))
print('------------------------------')
print('Simulation Completed')
print('------------------------------')
print('------------------------------')

######################################################################
# postprocessing #####################################################
t = np.array(ts)
sol = np.array(sol)
n = sol[:, 0]
c1 = sol[:, 1]
c2 = sol[:, 2]
c3 = sol[:, 3]
c4 = sol[:, 4]
c5 = sol[:, 5]
c6 = sol[:, 6]

plt.figure(1)
plt.plot(t, n)
plt.xlabel("Time [s]")
plt.ylabel("c1 [L^-3]")
plt.title("Neutron concentration")
plt.grid()

plt.figure(2)
plt.plot(t, c1)
plt.xlabel("Time [s]")
plt.ylabel("c1 [L^-3]")
plt.title("First group of precursors concentration")
plt.grid()

plt.figure(3)
plt.plot(t, c2)
plt.xlabel("Time [s]")
plt.ylabel("c2 [L^-3]")
plt.title("Second group of precursors concentration")
plt.grid()

# to show all figures
#plt.show()
#