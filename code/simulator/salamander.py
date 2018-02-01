# main page for salamander
from bot import *
import scipy.io
# global variables
# dutyf, interPhi, meterPhi, outerPhi


def main():
    Aleg = math.pi/6
    beta0 = np.array([1, -1, 1, -1])*math.pi/2
    leg_act = lambda x: np.array([leg_act_col(t, dutyf), 
                                leg_act_col(t+interPhi, dutyf),
                                leg_act_col(t+meterPhi, dutyf), 
                                leg_act_col(t+outerPhi, dutyf)]).flatten()
    beta = lambda x: np.array([leg_angle_col(Aleg, x, dutyf), 
                            leg_angle_col(Aleg, x+interPhi, dutyf),
                            leg_angle_col(Aleg, x+meterPhi, dutyf),
                            leg_angle_col(Aleg, x+outerPhi, dutyf)]).flatten()\
                            + beta0
    beta_0 = lambda x: np.array([-leg_angle_col(Aleg, x, dutyf), 
                            leg_angle_col(Aleg, x+interPhi, dutyf),
                            -leg_angle_col(Aleg, x+meterPhi, dutyf),
                            leg_angle_col(Aleg, x+outerPhi, dutyf)]).flatten()
    t_ini = t = 0
    T = 100
    dt = 2*math.pi/T
    f_num = 0
    alpha_pre_col, leg_col, beta_col = [], [], []
    xi = np.array([-math.pi/2, 0, 0, l_d+l_e/2, l_d+l_e/2, 
                l_d+l_e/2, l_d+l_e/2])
    data = scipy.io.loadmat('g_leg_initial.mat')
    g_leg_init = np.array(data['g_leg'])
    for i in range(5):
        while(t<t_ini + 4*math.pi):
            g_leg, xi, slip = get_config(g_leg, beta(t), alpha(t))






if __name__ == '__main__':
    main()