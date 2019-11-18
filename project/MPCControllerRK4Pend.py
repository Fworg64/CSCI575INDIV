import numpy as np
import scipy.optimize

from . import RK4Pendulum

class MCPControllerRK4Pend:
    def __init__(theta0, omega0, n_steps = 3, control_dt=.2, sim_dt=.02):
        self.states = [theta0, omega0]
        self.torque_plan = [0] * n_steps
        self.n_steps = n_steps
        self.control_dt = control_dt
        self.sim_dt = sim_dt

    def mpc_get_cost(self, pendulum_states, des_theta, torque_plan, cost_params=None):
        pend = RK4Pendulum(pendulum_states[0], pendulum_states[1], 
                           pendulum_states[2], sim_dt)
        frame_time = self.n_steps * self.control_dt
        sim_steps = int(self.control_dt / self.sim_dt)
        theta_rec = [0] * sim_steps * self.n_steps
        omega_rec = [0] * sim_steps * self.n_steps
        alpha_rec = [0] * sim_steps * self.n_steps
        for index in range(0, self.n_steps):
            pend.set_input_torque(torque_plan[index])
            [theta_rec[(index-1)*sim_steps:index*sim_steps],
             omega_rec[(index-1)*sim_steps:index*sim_steps],
             alpha_rec[(index-1)*sim_steps:index*sim_steps]] =  pend.simulate(self.control_dt)

        # Get cost vars ready
        des_theta_vec = [des_theta] * len(theta_rec)
        theta_arr = np.array(theta_rec)
        des_theta_arr = np.array(des_theta_vec)
        cost_arr = np.array(theta_arr - des_theta_arr)
        if cost_params=None:
            weight_arr = np.ones((len(theta_rec),))
        weight_arr = np.diag(weight_arr)
        # Do cost calculation with params
        precost = weight_arr @ cost_arr
        cost = np.transpose(cost_arr) @ precost
        return cost


    def mpc_force(self, pendulum_states, des_theta, cost_params, initial_torque_plan=None):
        # Minimize mpc get cost with given cost_params
      if initial_torque_plan is None:
          self.torque_plan = [0] * self.n_steps
      else:
          self.torque_plan = initial_torque_plan
      
      def obj(torque_plan):
          return self.mpc_get_cost(pendulum_states, des_theta, torque_plan, cost_params)
      torque_plan_out = scipy.optimize.minimize(obj, self.torque_plan, method = 'CG')
      print(torque_plan_out)
      self.torque_plan = torque_plan_out
      return torque_plan_out
      
