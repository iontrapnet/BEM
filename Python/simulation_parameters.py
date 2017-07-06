from __future__ import division
import numpy as np
from MField import MFieldInit

class simulation_parameters(object):

    def __init__(self):

        self.number_ions = 1
        # trap frequencies
        self.f_drive = 30.0 * 10**6  # Hz
        self.f_x = 4.0 * 10**6  # Hz
        self.f_y = 3.0 * 10**6  # Hzfield_type
        self.f_z = 0.2 * 10**6  # Hz
        self.field_type = 1
        if self.field_type == 1:
            MFieldInit(
                '..\\Model\\4rod\\167634622912717531', 
                [-0.005,0.005,100], 
                [-0.005,0.005,100], 
                [2.095,2.105,100], 
                '5d49bb9d6704e98ea598bb82a952b646'
                )
            self.f_drive = 12.0 * 10**6
            self.u_ac = 500
            self.u_dc = 30
            self.q_coeff = 1e2 * 1.6021766e-19
        # simulation parameter
        self.damping = 0  # optional velocity damping, useful for finding equlibrium positions
        self.simulation_duration = 0.002  # seconds
        self.timestep = (1 / self.f_drive) / 100  # seconds
        # ion parameters
        self.atomic_unit = 1.6605402e-27
        self.mass = 40 * 1.6605402e-27  # 40 amu in kg
        # k =  U.e**2 / (4.0 * U.pi * U.eps0)
        self.coulomb_coeff = 2.30707955552e-28
        self.hbar = 1.05457266913e-34
        self.transition_gamma = (1 / (7.1 * 10**-9))  # Gamma = 1 / Tau
        self.wavelength = 397 * 10**-9
        self.transition_k_mag = 2 * np.pi / (397 * 10**-9)
        # cooling laser
        self.cooling_on = False
        self.cooling_saturation = 0
        self.cooling_laser_detuning = -.5 * self.transition_gamma
        self.cooling_laser_direction = np.array([1., 1., 1.])
        self.cooling_laser_direction = self.cooling_laser_direction / \
            np.sqrt(np.sum(self.cooling_laser_direction**2))  # normalized
        self.cooling_laser_center = np.array([0.0, 0.0, 0.0])
        self.cooling_laser_waist = 0.001  # meters
        self.cooling_laser_pulsed = False
        # heating laser
        self.heating_on = False
        self.heating_saturation = 0
        self.heating_laser_detuning = .5 * self.transition_gamma
        self.heating_laser_direction = np.array([1., 1., 1.])
        self.heating_laser_direction = self.heating_laser_direction / \
            np.sqrt(np.sum(self.heating_laser_direction**2))  # normalized
        self.heating_laser_center = np.array([0.0, 0.0, 0.0])
        self.heating_laser_waist = 0.001  # meters
        self.heating_laser_pulsed = False

    @property
    def total_steps(self):
        return int(self.simulation_duration / self.timestep)
