# main page for salamander
from bot import *
import scipy.io
from math import sin, cos, pi
# global variables
# dutyf, interPhi, meterPhi, outerPhi

###########################################################
#           Simulator Models
###########################################################
def alpha_func(params, t):
    '''
    Body angles
    :return: np.array(size = (1,3), dtype = float)
    '''
    alpha = np.array([params.a2_1 * sin(t + params.b2_1),
                      params.a2_2 * sin(t + params.b2_2),
                      0.0])
    return alpha

def leg_act_func(params, t):
    '''
    Leg activation function
    :return: np.array(size = (1,4), dtype = float)
    '''
    leg_act = np.array([f_leg_act(t, params.dutyf)[0],
                        f_leg_act(t + params.interPhi, params.dutyf)[0],
                        f_leg_act(t + params.meterPhi, params.dutyf)[0],
                        f_leg_act(t + params.outerPhi, params.dutyf)[0]],
                        dtype = int)
    return leg_act

def beta_func(params, t):
    '''
    Leg movement
    :return:
    '''
    beta0 = np.array([1, -1, 1, -1], dtype = np.float) * (pi / 2)
    beta1 = np.array([ params.c_FR * f_leg(params.Aleg, t, params.dutyf)[0],
                     - params.c_FL * f_leg(params.Aleg, t + params.interPhi, params.dutyf)[0],
                       params.c_HR * f_leg(params.Aleg, t + params.meterPhi, params.dutyf)[0],
                     - params.c_HL * f_leg(params.Aleg, t + params.outerPhi, params.dutyf)[0]],
                     dtype = np.float)
    beta = beta0 + beta1
    return beta

def salamander_move(params):
    '''
    Run the simulation and return the displacement value
    :return: displacement
    '''
    # Initial positions for all legs (g_leg_initial.mat)
    data = scipy.io.loadmat('g_leg_initial.mat')
    g_leg = data['g_leg'][0] # Shape = (4,)

    # Initial head position
    head_pos_init = np.array([0.0, 0.0], dtype = np.float)

    # Get leg constants
    (l_a, l_b, l_c, l_d, l_e) = params.length_constants

    # Initial body frame position
    x_i = np.array([- pi/2, 0.0, 0.0, l_d+l_e/2.0, l_d+l_e/2.0, l_d+l_e/2.0, l_d+l_e/2.0])

    # Move the salamander
    t = params.t0

    while (t < params.t0 + 2 * pi):
        alpha = alpha_func(params, t).tolist()
        beta = beta_func(params, t).tolist()
        leg_act = leg_act_func(params, t).tolist()

        # print('\n************** I = {}, leg_act = {}'.format(i, leg_act))

        g_leg, x_i, slip = get_config(g_leg, beta, alpha, x_i, leg_act, params.length_constants)

        # Wrap up angle change x_i[0]
        if (x_i[0] < - pi):
            x_i[0] += pi
        elif (x_i[0] > pi):
            x_i[0] -= pi

        # print('\n************** After: ****************')
        # print("g_leg = \n{}".format(g_leg))
        # print("x_i = {}".format(x_i.tolist()))

        t += params.dt


    # Calculate Results
    head_pos_end = np.array([x_i[1], x_i[2]])
    displacement = np.linalg.norm(head_pos_end - head_pos_init)
    angleChange = x_i[0] + pi/2


    return (displacement, angleChange)


###########################################################
#           Simulator Driver
###########################################################
def init_constants(params):
    # Amplitude for legs
    params.Aleg = math.pi/3.0

    # Robot constants
    (l_a, l_b, l_c, l_d, l_e) = (2.25, 12.1, 4.55, 6.8, 10.35)
    params.length_constants = (l_a, l_b, l_c, l_d, l_e)

    # Time constants
    T = 100 # Number of time steps
    params.dt = 2 * pi / T

    return params

def run_baximander(dutyf = 0.5,
                   interPhi = pi, meterPhi = pi, outerPhi = pi,
                   a2_1 = 0.1, a2_2 = 0.1, b2_1 = 0.1, b2_2 = 0.1,
                   c_FL = 0.1, c_FR = 0.1, c_HL = 0.1, c_HR = 0.1):

    class Struct(object): pass
    params = Struct()

    # Set paramsters
    params.dutyf = dutyf
    params.interPhi = interPhi
    params.meterPhi = meterPhi
    params.outerPhi = outerPhi
    params.a2_1 = a2_1
    params.a2_2 = a2_2
    params.b2_1 = b2_1
    params.b2_2 = b2_2
    params.c_FL = c_FL
    params.c_FR = c_FR
    params.c_HL = c_HL
    params.c_HR = c_HR

    # Other parameters
    params.t0 = 0.0
    params = init_constants(params)

    (displacement, angleChange) = salamander_move(params)

    return (displacement, angleChange)

def test():
    dutyf = 0.9
    interPhi = 4.01212602265
    meterPhi = 1.18943521037
    outerPhi = -1.41714577328
    a2_1 = 0.701045249749
    a2_2 = 0.905461661386
    b2_1 = 2.86276125191
    b2_2 = 0.0474152993261
    c_FL = 0.900959302491
    c_FR = 0.790483551065
    c_HL = 0.84577315574
    c_HR = 0.606085753847

    (displacement, angleChange) = run_baximander(dutyf, interPhi, meterPhi, outerPhi,
                                                 a2_1, a2_2, b2_1, b2_2, c_FL, c_FR, c_HL, c_HR)

    print('Displacement = {} cm'.format(displacement))
    print('Angle change = {} rad'.format(angleChange))

if __name__ == '__main__':
    test()




