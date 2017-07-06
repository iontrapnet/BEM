'''
The function is designed to invoke the python program to simulate multiple ion dynamics with idealized quadropole potential, Coulomb interaction between ions, laser cooling and laser heating
'''
from simulation import simulator
import numpy as np
from simulation_parameters import simulation_parameters
from equilbrium_positions import equilibrium_positions as equil
from matplotlib import pyplot

# import os
# os.chdir(os.path.dirname(__file__))

def dynamics_simulation(configuration):
    '''
    The function requires a parameter configuration, which is a list or numpy array

        configuration = [trap_configuration, operation1, operation2, operation3...]

    trap_configuration consists of frequencies setting, ion settings

        trap_configuration = [number_ions, frequency_drive, frequency_x, frequency_y, frequency_z, damping, mass, transition_gamma, wavelength, use_harmonic_approximation, starting_positions, starting_positions]

        e.g. trap_configuration = [1, 30.0*10**6, 4.0 * 10**6, 3.0 * 10**6, 0.2 * 10**6, 0, 40, 1 / (7.1 * 10**-9), 397*10**-9, False, [[0,0, z1],[0,0,z2]], [[0,0,0],[0,0,0]]]

        z1, z2 are equilibrium positions obtained by equilbrium_positions.py, and you can assign a configuration without starting_positions and starting_velocity, which can be generated by the program, see below as an example.

        trap_configuration = [2, 30.0*10**6, 4.0 * 10**6, 3.0 * 10**6, 0.2 * 10**6, 0, 40, 1 / (7.1 * 10**-9), False]

        In details, the example is actually the default settings for the simulation, testing 40Ca.

        You should note that number_ions should be an integer, frequencies should all be double in Hz, damping is a dimensionless quantity, transition_gamma is 1/tau, use_harmonic_approximation should be assigned True if you want to deal with the simulation in pseudo potential approximation

        If you just want to run this simulation with those default settings, please assign trap_configuration = 0

        You can also change the number of default ion with a non-zero integer to tran_configuration

    operation indicates a duration with some laser interactions or not

        operation = [simulation_duration, cooling_laser, heating_laser,iteration_times, timesteps]

        laser = [saturation, laser_detuning, laser_direction, laser_center, laser_waist, laser_pulsed]

        In operation, four parameters are required
            simulation_duration = (double) time, unit in seconds
            iteration_times = (int) n, it is used to run the same simulation multiple times and average the results

        cooling_laser and heating_laser share the same data structure as laser above.

            saturation = 2Omega^2/Gamma^2

            laser_detuning = (double) detuning, is a factor compared to Gamma, for example if you want to set laser_detuning = 0.5 * Gamma, you actually just need to set laser_detuning = -0.5, but do notice that laser_detuning can be positive for heating, and you must not assign wrong values to lasers, cooling requires detuning to be negative, but positive for heating

            laser_direction = [x,y,z], not necessary to make it normalized

            laser_center = [x,y,z], a point in space

            laser_waist = (double) waist, in meters

            laser_pulsed = True/False, indicating pulsed laser or not
        **************************************************************************
        P.S. if any kind of laser is not applied, please assign laser to be 0, heating_laser = 0 for example, operation = [simulation_duration, iteration_times, cooling_laser, 0]
        !!! Numpy array is prefered for array
        **************************************************************************

        timesteps is not necessary, the default would be 100, which means timestep in the parameters 1/f_drive/100, but you can try it for some different numbers,
        for example, timesteps = 5, cause simulation timestep = 1/f_drive/5

    If you want to applay different laser or anything, please add more operations, but the least number of operations is 1, even you don't apply any laser, just testing its dynamics caused by the quadropole potential and Coulomb interation.

    Returns a list consisting [positions, excitations, time_steps]

        positions = np.array([positions_1,positions_2......])
        excitations = np.array([excitations_1,excitations_2......])
        time_steps = np.array([timesteps_1,timesteps_2......])
    '''
    parameters = simulation_parameters()
    if configuration[0] == 0:
        starting_positions = np.zeros((parameters.number_ions, 3))
        starting_positions[:, 0] = parameters.number_ions * [0]  # x
        starting_positions[:, 1] = parameters.number_ions * [0]  # y
        starting_positions[:, 2] = equil.get_positions(
            parameters.number_ions, parameters.f_z, parameters)
        starting_positions += (np.random.random(
            (parameters.number_ions, 3)) - .5) * 1e-6
        starting_velocities = np.zeros((parameters.number_ions, 3))
    elif type(configuration[0]) == type(1) and configuration[0] > 0:
        parameters.number_ions = configuration[0]
        starting_positions = np.zeros((parameters.number_ions, 3))
        starting_positions[:, 0] = parameters.number_ions * [0]  # x
        starting_positions[:, 1] = parameters.number_ions * [0]  # y
        starting_positions[:, 2] = equil.get_positions(
            parameters.number_ions, parameters.f_z, parameters)
        starting_positions += (np.random.random(
            (parameters.number_ions, 3)) - .5) * 1e-6
        starting_velocities = np.zeros((parameters.number_ions, 3))
    elif len(configuration[0]) == 10:
        parameters.number_ions = configuration[0][0]
        parameters.f_drive = configuration[0][1]
        parameters.f_x = configuration[0][2]
        parameters.f_y = configuration[0][3]
        parameters.f_z = configuration[0][4]
        parameters.damping = configuration[0][5]
        parameters.mass = configuration[0][6] * parameters.atomic_unit
        parameters.transition_gamma = configuration[0][7]
        parameters.wavelength = configuration[0][8]
        parameters.transition_k_mag = 2 * np.pi / parameters.wavelength
        parameters.use_harmonic_approximation = configuration[0][9]
        starting_positions = np.zeros((parameters.number_ions, 3))
        starting_positions[:, 0] = parameters.number_ions * [0]  # x
        starting_positions[:, 1] = parameters.number_ions * [0]  # y
        starting_positions[:, 2] = equil.get_positions(
            parameters.number_ions, parameters.f_z, parameters)
        starting_positions += (np.random.random(
            (parameters.number_ions, 3)) - .5) * 1e-6
        starting_velocities = np.zeros((parameters.number_ions, 3))
    elif len(configuration[0]) == 12:
        parameters.number_ions = configuration[0][0]
        parameters.f_drive = configuration[0][1]
        parameters.f_x = configuration[0][2]
        parameters.f_y = configuration[0][3]
        parameters.f_z = configuration[0][4]
        parameters.damping = configuration[0][5]
        parameters.mass = configuration[0][6] * parameters.atomic_unit
        parameters.transition_gamma = configuration[0][7]
        parameters.wavelength = configuration[0][8]
        parameters.transition_k_mag = 2 * np.pi / parameters.wavelength
        parameters.use_harmonic_approximation = configuration[0][9]
        starting_positions = np.array(configuration[0][10])
        starting_velocities = np.array(configuration[0][11])
    else:
        raise Exception("Wrong input for trap settings")
    starting_time = 0
    starting_excitations = np.zeros(parameters.number_ions, dtype=np.uint8)
    excitations = []
    positions = []
    time = []
    order = 0
    for operation in configuration[1:]:
        order += 1
        print "!!!Operation", order, "starts!"
        parameters.simulation_duration = operation[0]
        iteration_times = 1
        if len(operation) == 5:
            iteration_times = operation[3]
            parameters.time_steps = (1 / parameters.f_drive) / operation[4]
        if len(operation) == 4:
            iteration_times = operation[3]
        chunksize = int(parameters.total_steps * 3e-3)
        # initialize cooling laser
        if not operation[1] == 0:
            parameters.cooling_on = True
            parameters.cooling_saturation = float(operation[1][0])
            parameters.cooling_laser_detuning = operation[
                1][1] * parameters.transition_gamma
            parameters.cooling_laser_direction = np.array(
                operation[1][2]).astype(float)
            parameters.cooling_laser_direction = parameters.cooling_laser_direction / \
                np.sqrt(np.sum(parameters.cooling_laser_direction**2)
                        )  # normalized
            parameters.cooling_laser_center = np.array(
                operation[1][3]).astype(float)
            parameters.cooling_laser_waist = float(operation[1][4])
            parameters.cooling_laser_pulsed = operation[1][5]
        # initialize heating laser
        if not operation[2] == 0:
            parameters.heating_on = True
            parameters.heating_saturation = float(operation[2][0])
            parameters.heating_laser_detuning = operation[
                2][1] * parameters.transition_gamma
            parameters.heating_laser_direction = np.array(
                operation[2][2]).astype(float)
            parameters.heating_laser_direction = parameters.heating_laser_direction / \
                np.sqrt(np.sum(parameters.heating_laser_direction**2)
                        )  # normalized
            parameters.heating_laser_center = np.array(
                operation[2][3]).astype(float)
            parameters.heating_laser_waist = float(operation[2][4])
            parameters.heating_laser_pulsed = operation[2][5]
        positions_i = np.zeros(
            (parameters.number_ions, parameters.total_steps, 3))
        excitations_i = np.zeros(
            (parameters.total_steps, parameters.number_ions), dtype=np.uint8)
        for i in range(iteration_times):
            print 'ITERATION', i + 1
            simulator_i = simulator(parameters)
            positions_operation_i, excitations_operation_i = simulator_i.simulation(
                starting_positions, starting_velocities, starting_excitations, random_seeding=i)
            positions_i += positions_operation_i
            excitations_i += excitations_operation_i
        print "!!!Operation", order, "ends"
        positions_i /= iteration_times
        excitations_i = np.round(excitations_i / iteration_times)
        timesteps = np.arange(parameters.total_steps) * \
            parameters.timestep * 10**6 + starting_time
        # update the starting time, positions and excitations for next
        # operation
        starting_time = timesteps[-1]
        starting_velocities = (
            positions_i[:, -1, :] - positions_i[:, -2, :]) / parameters.timestep
        starting_positions = positions_i[:, -1, :]
        starting_excitations = excitations_i[-1, :]
        positions.append(positions_i)
        excitations.append(excitations_i)
        time.append(timesteps)
        positions_chunk = positions_i[:, :(parameters.total_steps // chunksize)
                                      * chunksize, :].reshape((parameters.number_ions, -1, chunksize, 3))
        timesteps_chunk = timesteps[
            :(parameters.total_steps // chunksize) * chunksize].reshape((-1, chunksize))
        positions_center = positions_chunk.mean(axis=2) * 10**6
        timesteps_center = timesteps_chunk.mean(axis=1)
        # positions_center = positions_i * 10**6
        # timesteps_center = timesteps
        for num in range(parameters.number_ions):
            pyplot.figure((order - 1) * parameters.number_ions + num + 1)
            pyplot.title("Trajectory of ion {} in operation {}".format(
                num + 1, order), fontsize=12)
            pyplot.xlabel(r'Time($\mu s$)', fontsize=10)
            pyplot.ylabel(r'Positions($\mu m$)', fontsize=10)
            pyplot.tick_params(axis='both', labelsize=6)
            pyplot.plot(timesteps_center, positions_center[
                        num, :, 0], 'blue', label='x')
            pyplot.plot(timesteps_center, positions_center[
                        num, :, 1], 'red', label='y')
            pyplot.plot(timesteps_center, positions_center[
                        num, :, 2], 'black'
, label='z')
            pyplot.grid(True, which='both')
            pyplot.legend()
            filename = "ion " + str(num + 1) + \
                " in operation " + str(order) + ".png"
            pyplot.savefig(filename, dpi=100)
    pyplot.show()
    positions = np.array(positions)
    excitations = np.array(excitations)
    time = np.array(time)
    return [positions, excitations, time]
if __name__ == '__main__':
    # Test 1 for the default settings without laser, please run only a single
    # case below, due to the figures maybe sharing same name, choose one and
    # make others commented with #

    # 1 ca+ ion without laser
    # parameters position, excitation, time can be used for further analysis
    # Defaultly, we generate some figures for every individual ion to show its
    # positions
    position, excitation, time = dynamics_simulation(
        [0, [100e-6, 0, 0]])

    # 1 ion with 2 operations, but we don't applay any laser, the second
    # operation runs two times for average values over two individual
    # simulations with the same initial conditions
    # position, excitation, time = dynamics_simulation(
    #     [0, [100e-6, 0, 0], [100e-6, 0, 0, 2]])

    # 5 ion with 1 operation and no laser, a non-zero value for the first
    # argument will change the number of ions, and 50 in the second argument
    # will change the timestep into 1/f_drive/50, also you can add more
    # operations
    # position, excitation, time = dynamics_simulation(
    #     [5, [100e-6, 0, 0, 1, 50]])

    # # non-trivial construction, assign the properties of trap manually
    # # 2 ion of 40 ca without laser, starting positions and starting velocities
    # # can be neglected
    # trap_configuration = [2, 30.0 * 10**6, 4.0 * 10**6, 3.0 * 10**6, 0.2 * 10**6, 0, 40,
    #                       1 / (7.1 * 10**-9), 397 * 10**-9, False, [[4e-6, 0, -5e-5], [0, 3e-5, 5e-5]], [[0, 0, 0], [0, 0, 0]]]
    # print len(trap_configuration)
    # operation1 = [100e-6, 0, 0, 2, 100]
    # operation2 = [100e-6, 0, 0, 3, 50]
    # position, excitation, time = dynamics_simulation(
    #     [trap_configuration, operation1, operation2])

    # # 5 ion of 40 ca with laser
    # trap_configuration = [5, 30.0 * 10**6, 4.0 * 10**6, 3.0 * 10**6, 0.2 * 10**6, 0, 40,
    #                       1 / (7.1 * 10**-9), 397 * 10**-9, False]
    # # operation 1 with no heating laser
    # operation1 = [100e-6, [3.0, -0.5, [1., 1., 1.],
    #                        [0, 0, 0], 0.001, False], 0]
    # # operation 2 with a pulsed heating laser, confirm detuning for cooling
    # # should be negative but positive for heating
    # operation2 = [100e-6, [5.0, -0.5, [3., 1., 1.], [5, 2, 1], 0.001,
    #                        False], [0.1, 0.5, [3., 1., 1.], [5, 2, 1], 0.001, True], 2, 100]
    # position, excitation, time = dynamics_simulation(
    #     [trap_configuration, operation1, operation2])
