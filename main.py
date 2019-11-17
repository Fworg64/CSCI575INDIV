# dTheta = w
# w+ = w + (l1 + 2*l2 + 2*l3 + l4)*dt / 6
# th+= th + (k1 + 2*k2 + 2*k3 + k4)*dt / 6
# l1 = g(t, w; th) = dw
# l2 = g(t + dt/2, wn + l1*dt/2)
# k1 = f(t, th; w) = dth

import matplotlib.pyplot
import numpy as np
import scipy.optimize

from datetime import datetime  

now = datetime.now()

cont_dt = .2
dt = .02
dt_2 = dt/2.0


def MPC_get_cost(theta0, omega0, des_theta, torque_plan, control_dt, sim_dt, weight_arr=None):
  # Calculates position - des_pos at each stage and takes RMS
  # sim_dt <10x< control_dt
  theta = theta0
  omega = omega0
  omega_rec = [omega,]
  theta_rec = [theta,]
  force_rec = [torque_plan[0],]
  indices = [0,]
  
  def force(t): # interpolate torque plan with control_dt and sim_dt
    # each tp entry is control_dt long
    steps = t/sim_dt
    ratio = control_dt/sim_dt
    index = np.min([int(np.floor((steps-1)/ratio)), len(torque_plan)-1])
    #print("ehad")
    #print(index)
    #print(len(torque_plan))
    return torque_plan[index]

  def dtheta(t, theta, omega): # parameterized dynamics
    value = omega # Just integrate it
    return value

  def domega(t, theta, omega): # parameterized dynamics
    value = -np.sin(theta) - .5*omega + force(t)
    return value

  for t in np.arange(0.0,10,dt):
    # RK4 Numerical Integration for simulation
    k1 = dtheta(t, theta, omega)
    l1 = domega(t, theta, omega)

    k2 = dtheta(t + dt_2, theta + k1*dt_2, omega + l1*dt_2)
    l2 = domega(t + dt_2, theta + k1*dt_2, omega + l1*dt_2)

    k3 = dtheta(t + dt_2, theta + k2*dt_2, omega + l2*dt_2)
    l3 = domega(t + dt_2, theta + k2*dt_2, omega + l2*dt_2)

    k4 = dtheta(t + dt, theta + k3*dt, omega + l3*dt)
    l4 = domega(t + dt, theta + k3*dt, omega + l3*dt)

    theta = theta + (k1 + 2*k2 + 2*k3 + k4)*dt/6
    omega = omega + (l1 + 2*l2 + 2*l3 + l4)*dt/6
    # record states
    omega_rec.append(omega)
    theta_rec.append(theta)
    force_rec.append(force(t))
    indices.append(indices[-1]+1)
  des_theta_vec = [des_theta] * len(theta_rec) 
  theta_arr = np.array(theta_rec)
  des_theta_arr = np.array(des_theta_vec)
  cost_arr = np.array(theta_arr - des_theta_arr)
  cost_arr = cost_arr.reshape((len(theta_rec),1))
  if weight_arr is None:
    weight_arr = np.ones((len(theta_rec),))
  #print(cost_arr.shape)
  weight_arr = np.diag(weight_arr)
  #print(weight_arr.shape)
  #print(weight_arr)
  precost = weight_arr @ cost_arr
  #print(precost.shape)
  #print(precost)

  cost = np.transpose(cost_arr) @ precost
  #double_cost = float(cost)
  #cost = sum([(theta - des_theta)**2 for theta in theta_rec])
  return cost[0]
  

def MPC_force(pos, vel, des_pos):
  # Set up variables
  torque_plan = [0,0,0]
  control_dt = cont_dt;
  sim_dt = dt;
  def obj(torque_plan):
    return MPC_get_cost(pos, vel, des_pos, torque_plan, control_dt, sim_dt)
    
  torque_plan = scipy.optimize.minimize(obj, torque_plan, method = 'CG')
  print(torque_plan)
  return torque_plan.x[0]
  
  # Populate Torque plan and get initial cost
  #Loop A (while cost > tol AND curr_iter < max_iter ):
  # Solve for current torque plan and evaluate cost
  
  # Adjust torque plan to minimize cost
  # End Loop A



def dtheta(t, theta, omega):
  value = omega # Just integrate it
  return value

def domega(t, theta, omega, control_force):
  value = -np.sin(theta) - .5*omega + control_force
  return value

theta = .2
omega = 0

omega_rec = [omega,]
theta_rec = [theta,]
force_rec = [0,]
indices = [0,]

frames = 0;

for t in np.arange(0.0,1.2,dt):
  
  # Periodically evaluate control
  if t >= frames*cont_dt:
    frames = frames+1
    control_force = MPC_force(theta, omega, .6)
  
  # RK4 Numerical Integration for simulation
  k1 = dtheta(t, theta, omega)
  l1 = domega(t, theta, omega, control_force)
  
  k2 = dtheta(t + dt_2, theta + k1*dt_2, omega + l1*dt_2)
  l2 = domega(t + dt_2, theta + k1*dt_2, omega + l1*dt_2, control_force)
  
  k3 = dtheta(t + dt_2, theta + k2*dt_2, omega + l2*dt_2)
  l3 = domega(t + dt_2, theta + k2*dt_2, omega + l2*dt_2, control_force)
  
  k4 = dtheta(t + dt, theta + k3*dt, omega + l3*dt)
  l4 = domega(t + dt, theta + k3*dt, omega + l3*dt, control_force)
  
  theta = theta + (k1 + 2*k2 + 2*k3 + k4)*dt/6
  omega = omega + (l1 + 2*l2 + 2*l3 + l4)*dt/6
  # record states
  omega_rec.append(omega)
  theta_rec.append(theta)
  force_rec.append(control_force)
  indices.append(indices[-1]+1)
  
# Get Datestamp
then = datetime.now()
time_take = then - now
print(str(time_take) + " seconds")
now_str = now.strftime("%m_%d_%Y_%H_%M_%S")
# Plot states
matplotlib.pyplot.ioff()

fig = matplotlib.pyplot.figure()
matplotlib.pyplot.plot(indices, theta_rec, indices, omega_rec)
matplotlib.pyplot.legend(["theta", "omega"])
matplotlib.pyplot.savefig('./figs/theta_omega_plot' + now_str + '.png')
matplotlib.pyplot.close(fig)

fig = matplotlib.pyplot.figure()
matplotlib.pyplot.plot(indices, force_rec)
matplotlib.pyplot.legend(["force"])
matplotlib.pyplot.savefig('./figs/force_plot' + now_str + '.png')
matplotlib.pyplot.close(fig)

print("Goodbye!")

# Change forcing equation to MPC solver which uses same model but smaller 
# Control time step
