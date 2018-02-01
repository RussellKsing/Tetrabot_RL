from bot import *
import scipy.io
import sys
from pprint import pprint
#leg function, *inhead
def test_1():
    print("Start Testing...")
    print("Output:")
    #leg_act_col(time, dutyf): why is time -1?
    print("Funtion leg_act_col: ", leg_act_col([1, -1, 10], 0.75))
    print("Funtion leg_angle_col: ", leg_angle_col(10, [1, -1, 10], 0.75))
    print("Matlab output: 9.8978e+00  -9.8978e+00  -3.7415e+00")
    print("Function R: ", R(0.5*math.pi))
    print("""Matlab output: 6.1232e-17,-1.0000e+00,0;1.0000e+00, 6.1232e-17, 0; 
                0, 0, 1.0000e+00""")
    # print("Function joints_in_head: ", joints_in_head(np.array([[1, -1, 10],[1, -1, 10]]), 10))
    pprint("Function joints_in_head: ")
    pprint (joints_in_head(np.array([1, -1, 10]), 10))
    alpha = np.array([1, -1, 10]);
    beta = np.array([2,-3, 4, 6]);
    hind = 2
    g, gFL, gFR, gHL, gHR = frames_in_head(alpha, beta, hind)
    # sys.exit()
    print("\nFunction frames_in_head: ")
    pprint("\ng: ", g)
    print("\n gFL: ", gFL)
    # gFL =

    #   -4.1615e-01  -9.0930e-01  -3.9066e+00
    #    9.0930e-01  -4.1615e-01   1.3086e+01
    #             0            0   1.0000e+00
    print("\ngFR: ", gFR)
      # -9.8999e-01   1.4112e-01  -9.2936e+00
      # -1.4112e-01  -9.8999e-01  -5.8748e+00
      #           0            0   1.0000e+00
    print("\ngHL: ", gHL)
      # -6.5364e-01   7.5680e-01  -6.1361e+00
      # -7.5680e-01  -6.5364e-01  -2.5545e+00
      #           0            0   1.0000e+00
    print("\ngHR: ", gHR)
    # gHR =

    #    9.6017e-01   2.7942e-01   9.0136e+00
    #   -2.7942e-01   9.6017e-01  -7.1730e+00
    #             0            0   1.0000e+00
    print(".....Section I Test Passed")


def test_2():
    print("Start Test Section II......")
    # test function G(a, x, y)
    # print(G(1,2,3))
    # passed
    xi = np.array([-1.5, -1.7, -9, 9, 9.3, 9.3, 9.3]);
    alpha = [1, -1, 10];
    beta = [2, -3, 4, 6];
    data = scipy.io.loadmat('g_leg_initial.mat')
    # sys.exit()
    g_leg = data['g_leg'][0]
    # print("g_leg:", g_leg)
    print(np.linalg.norm(g_leg[3][0:2,2]- g_leg[0][0:2,2]))
    t = 0
    lfs=0.25;
    interPhi=math.pi;
    meterPhi=lfs*2*math.pi;
    outerPhi=math.pi+meterPhi;
    leg_act_0 = [0, 1, 1, 1]
    g_leg, xi, slip = get_config(g_leg, beta, alpha, xi, leg_act_0)
    print("g_leg: ", g_leg)
    print("xi: ", xi)
    print("slip: ", slip)
    data_mat = scipy.io.loadmat('test_data.mat')
    # print(data_mat)
    g_leg_m, xi_m, slip_m = data_mat["g_leg"], data_mat["xi"], data_mat["slip"]
    print("g_leg_matlib: ", g_leg_m)
    print("xi_matlib: ", xi_m, "Function value: 5.511483895407483")
    print("slip_matlib: ", slip_m)
    print("g_leg_differences 0 :", g_leg_m[0][0] - g_leg[0])
    print("g_leg_differences 1:", g_leg_m[0][1] - g_leg[1])
    print("g_leg_differences 2:", g_leg_m[0][2] - g_leg[2])
    print("g_leg_differences 3:", g_leg_m[0][3] - g_leg[3])
    print(".....Section II Test Underreview!")

def main():
    # test_1() 
    test_2()
    return 42

if __name__ == '__main__':
    main()



