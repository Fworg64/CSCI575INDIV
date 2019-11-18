import numpy as np


class RK4Pendulum:
    def __init__(self, theta0, omega0, torque0, sim_dt):
      self.states = [theta0, omega0, torque0]
      self._dt = sim_dt
      self._dt_2   = sim_dt / 2.0

    def set_input_torque(self, torque):
        self.states[2] = torque

    def _torque(self):
        return self.states[2]

    def _dtheta(self, t, theta, omega):
        return omega

    def _domega(self, t, theta, omega):
        return -np.sin(theta) - 0.5*omega + self._torque()

    def simulate(self, time):
      theta_rec = [self.states[0],]
      omega_rec = [self.states[1],]
      alpha_rec = [self._domega(0, self.states[0], self.states[1]),]

      for t in np.arange(0.0,time,self._dt):
      # RK4 Numerical Integration for simulation
        k1 = self._dtheta(t, self.states[0], self.states[1])
        l1 = self._domega(t, self.states[0], self.states[1])

        k2 = self._dtheta(t + self._dt_2, self.states[0] + k1*self._dt_2,
                                          self.states[1] + l1*self._dt_2)
        l2 = self._domega(t + self._dt_2, self.states[0] + k1*self._dt_2,
                                          self.states[1] + l1*self._dt_2)

        k3 = self._dtheta(t + self._dt_2, self.states[0] + k2*self._dt_2,
                                          self.states[1] + l2*self._dt_2)
        l3 = self._domega(t + self._dt_2, self.states[0] + k2*self._dt_2,
                                          self.states[1] + l2*self._dt_2)

        k4 = self._dtheta(t + self._dt, self.states[0] + k3*self._dt,
                                        self.states[1] + l3*self._dt)
        l4 = self._domega(t + self._dt, self.states[0] + k3*self._dt,
                                        self.states[1] + l3*self._dt)

        self.states[0] = self.states[0] + (k1 + 2*k2 + 2*k3 + k4)*self._dt/6
        self.states[1] = self.states[1] + (l1 + 2*l2 + 2*l3 + l4)*self._dt/6

        theta_rec.append(self.states[0])
        omega_rec.append(self.states[1])
        alpha_rec.append(self._domega(t+self._dt, self.states[0], self.states[1]))

      return [theta_rec, omega_rec, alpha_rec]




