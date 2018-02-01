import numpy as np
import math
from scipy.optimize import minimize, fmin_slsqp
import copy
import sys
################################# Section I ####################################
# for an angle print it in radian as well as degree
def print_angle(ang_rad):
    print (ang_rad, math.degrees(ang_rad))

# F_leg_act.m
# bring the radian between 0 and 2pi by adding or subtracting
# 2pi repeatedly
def radian_transfer(t):
    # print ('radian transfer')
    # print_angle (t)
    while(t < 0 or t > 2*math.pi):
        # why is it subtracting pi and not 2pi?
        t = t + 2*math.pi if t < 0 else t - 2*math.pi
    # print_angle (t)
    return t

# time: 1d array
# dutyf: scalar
# F_leg_act.m
# it returns an array containing an indicator if the corresponding "time"
# values are less than pi*(1+d). why?
# probably checks collision?
def leg_act_col(time, dutyf):
    # print ('time: ', time)
    # print ('dutyf: ', dutyf)
    #set ans to 0's of length of time
    try:
        ans = np.zeros(len(time))
    except:
        time = np.array([time])
        ans = np.zeros(len(time))
    # sets the answer array to some indicator function 
    # - need to look into theory of it
    for i in range(len(time)):
        t = radian_transfer(time[i])
        ans[i] = 1 if t < 2*math.pi - (1 - dutyf)*math.pi else 0
    return ans

# F_Leg
# probably checks colision?
# checks a few if-else conitions based on time
# and returns a function of aleg and duty factor
# aleg == amplitude of the leg?
def leg_angle_col(Aleg, time, dutyf):
    # create a variable ans and set it to 0's of length of the time variable
    try:
        ans = np.zeros(len(time))
    except:
        time = np.array([time])
        ans = np.zeros(len(time))
    #
    for i in range(len(time)):
        # change the range of the angle between 0 to 2pi
        t = radian_transfer(time[i])
        if t < (1 - dutyf)*math.pi:
            # print ('1st cond')
            leg_angle = Aleg*math.sin(t/(2*(1-dutyf)))
        elif t < (2*math.pi - (1 - dutyf)*math.pi):
            # print ('2nd cond')
            leg_angle = Aleg*math.cos((t-(1-dutyf)*math.pi)*(math.pi/(2*math.pi\
                                                        -2*(1-dutyf)*math.pi)))
        else:
            # print ('3rd cond')
            leg_angle = -Aleg*math.cos((t-2*math.pi+\
                                (1- dutyf)*math.pi)/2/(1-dutyf))
        ans[i] = leg_angle
    # print (ans)
    return ans


# alpha: array
# used in jointsInHead
# solves the equation for x in aT xT = bT where T means transpose
def find_xa_b(a, b):
    # print (np.transpose(a))
    # print (np.transpose(b))
    # print (np.linalg.solve(np.transpose(a), np.transpose(b)))
    # # np.linalg.solve here solves for x in aTxT = bT
    return np.transpose(np.linalg.solve(np.transpose(a), np.transpose(b)))

# Some computations regated to the matrices and "F" matrices
def joints_in_head(alpha, L):
    try:
        n = alpha.shape[1] + 1
    except:
        n = 2
    # g is a list of matrices?
    print (n)
    # sys.exit()
    g = []
    g.append(np.eye(3))
    # returns a matrix with (0,2) indexed value as L
    f = F(L)
    # g[0] now has the inverse of f
    g[0] = find_xa_b(f, g[0])
    # print (g[0])
    for i in range(1, n-1, 1):
        # print('in the loop')
        temp = np.matmul(g[-1], f)
        temp = np.matmul(temp, f)
        # potential bug here!! what if alpha is a matrix?
        # print (alpha[i-1].shape)
        temp = np.matmul(temp, R(alpha[i-1]))
        g.append(temp)
    # print (g)
    # sys.exit()
    # the above loop has multiplication twice with f, but here only once
    # indexing problem for alpha
    temp = np.matmul(g[-1], f)
    temp = np.matmul(temp, R(alpha[n-2]))
    g.append(temp)
    f = F(L/3)
    temp = (np.matmul(g[-1], f))
    temp = np.matmul(temp, f)
    g.append(temp)
    return g

# FramesInHead
def matrix_multiply(*args):
    ans = args[0]
    for i in range(1, len(args), 1):
        ans = np.matmul(ans, args[i])
    return ans

# def frames_in_head(alpha, beta, hind, L, Lleg):
#     try: 
#         n = alpha.shape[1] + 1
#     except:
#         n = 2
#     g = []
#     g.append(np.eye(3))
#     temp = np.matmul(np.matmul(g[-1], F(3.6)), R(alpha[0]))
#     g.append(np.matmul(temp, F(L)))
#     g.append(np.matmul(g[-1], F(2*L)))
#     g.append(np.matmul(np.matmul(g[-1], F(L)), R(alpha(2)), F(3.6)))
#     gFL = matrix_multiply(g[0], R(math.pi/2), F(3.1), R(-math.pi/2), R(beta[0]),
#                         F(Lleg))
#     gFR = matrix_multiply(g[0], R(-math.pi/2), F(3.1), R(math.pi/2), R(beta[1]),
#                         F(Lleg))
#     gHL = matrix_multiply(g(hind-2), R(math.pi/2), F(3.1), R(-math.pi/2), 
#                         R(beta[2]), F(LLeg))
#     gHL = matrix_multiply(g(hind-2), R(-math.pi/2), F(3.1), R(math.pi/2), 
#                         R(beta[3]), F(LLeg))
#     return g, gFL, gFR, gHL, gHR

def frames_in_head(alpha, beta, hind):
    try: 
        n = alpha.shape[0] + 1
    except:
        n = 2
    # load data - what are these values?
    l_a, l_b, l_c, l_d, l_e = 2.25, 12.1, 4.55, 6.8, 10.35
    g = []
    g.append(np.eye(3))
    g.append(matrix_multiply(g[-1], F(l_a), R(alpha[0]), F(l_b/4)))
    g.append(matrix_multiply(g[-1], F(l_b/2)))
    g.append(matrix_multiply(g[-1], F(l_b/4), R(alpha[1]), F(l_a)))
    gFL = matrix_multiply(g[0], R(math.pi/2), F(l_c), R(-math.pi/2), R(beta[0]),
                        F(l_d+l_e/4))
    gFR = matrix_multiply(g[0], R(-math.pi/2), F(l_c), R(math.pi/2), R(beta[1]),
                        F(l_d+l_e/4))
    gHL = matrix_multiply(g[hind-2], R(math.pi/2), F(l_c), R(-math.pi/2), 
                        R(beta[2]), F(l_d+l_e/4))
    gHR = matrix_multiply(g[hind-2], R(-math.pi/2), F(l_c), R(math.pi/2), 
                        R(beta[3]), F(l_d+l_e/4))
    return g, gFL, gFR, gHL, gHR

# rotation matrix
# a = angle to rotate by
def R(a):
    ans = np.array([[math.cos(a), -math.sin(a), 0],
           [math.sin(a), math.cos(a), 0],
           [0, 0, 1]])
    return ans 

# sets (0,2) indexed value in the matrix as L and returns the matrix
# forward matrix
def F(L):
    ans = np.array([[1, 0, L], [0, 1, 0], [0, 0, 1]])
    return ans

################################# Section II ###################################

def G(a, x, y):
    ans = np.array([[math.cos(a), -math.sin(a), x],
         [math.sin(a), math.cos(a), y], [0, 0, 1]])
    return ans

def rms(x):
    return np.sqrt(np.mean(np.square(x)))

# what are getf, getg
def get_grad(x, getf, getg):
    return getf(x), getg(x)

# g_leg: list
# beta: 
# alpha: 
# xi: [-1.5 1 1 1.2 1.2 1.2]
# leg_act: [1 1 1 1]
def get_config(g_leg,beta,alpha,xi,leg_act,colorSpace = 1):
    # load data
    l_a, l_b, l_c, l_d, l_e = 2.25, 12.1, 4.55, 6.8, 10.35
    xi_pre, slip = copy.copy(xi), 0
    FR, FL, HR, HL = g_leg[0], g_leg[1], g_leg[2], g_leg[3]
    air_leg_length = l_d + l_e/4
    if leg_act[0] == 1 & leg_act[3] == 1:
        flag = 1
        if leg_act[1] == 1: flag0 = 1
        if leg_act[2] == 1: flag0 = 2
    elif leg_act[1] == 1 & leg_act[2] == 1:
        flag = 2
        if leg_act[0] == 1: flag0 = 0
        if leg_act[3] == 1: flag0 = 3
    xi = func(xi_pre, l_a, l_b, l_c, l_d, l_e, alpha, beta, FR, FL, HL, HR, 
            flag)
    g_leg, xi = cond(g_leg, xi,l_a, l_b, l_c, air_leg_length, alpha, beta, flag)
    slip = cond0(g_leg, slip, FR, FL, HR, HL, flag0)
    return g_leg, xi, slip

def cond(g_leg, xi,l_a, l_b, l_c, air_leg_length, alpha, beta, flag):
    if flag == 1:
        g_leg[0] = matrix_multiply(G(xi[0], xi[1], xi[2]), R(math.pi/2), F(l_c),
                                R(-math.pi/2), R(beta[0]/2), F(xi[3]))
        g_leg[1] = matrix_multiply(G(xi[0],xi[1],xi[2]), R(-math.pi/2), F(l_c),
                                R(math.pi/2), R(beta[1]), F(air_leg_length))
        g_leg[2] = matrix_multiply(G(xi[0],xi[1],xi[2]), F(l_a), R(alpha[0]), 
                                F(l_b), R(alpha[1]), F(l_a), R(math.pi/2), 
                                F(l_c), R(-math.pi/2), R(beta[2]), 
                                F(air_leg_length))
        g_leg[3] = matrix_multiply(G(xi[0],xi[1],xi[2]), F(l_a), R(alpha[0]), 
                                F(l_b), R(alpha[1]), F(l_a), R(-math.pi/2), 
                                F(l_c), R(math.pi/2), R(beta[3]), F(xi[3]))
        xi = np.append(xi[0:4], [air_leg_length, air_leg_length, xi[3]])
    elif flag == 2:
        g_leg[0] = matrix_multiply(G(xi[0], xi[1], xi[2]), R(math.pi/2), F(l_c),
                                R(-math.pi/2), R(beta[0]), F(air_leg_length))
        g_leg[1] = matrix_multiply(G(xi[0],xi[1],xi[2]), R(-math.pi/2), F(l_c),
                                R(math.pi/2), R(beta[1]), F(xi[3]))
        g_leg[2] = matrix_multiply(G(xi[0],xi[1],xi[2]), F(l_a), R(alpha[0]), 
                                F(l_b), R(alpha[1]), F(l_a), R(math.pi/2), 
                                F(l_c), R(-math.pi/2), R(beta[2]), 
                                F(xi[3]))
        g_leg[3] = matrix_multiply(G(xi[0],xi[1],xi[2]), F(l_a), R(alpha[0]), 
                                F(l_b), R(alpha[1]), F(l_a), R(-math.pi/2), 
                                F(l_c), R(math.pi/2), R(beta[3]),
                                F(air_leg_length))
        xi = np.append(xi[0:3], [air_leg_length, xi[3], xi[3], air_leg_length])
    return g_leg, xi

def cond0(g_leg, slip, FR, FL, HR, HL, flag0):
    if flag0 == 0: slip += np.linalg.norm(g_leg[0][0:2, 2]-FR[0:2, 2])
    elif flag0 == 1: slip += np.linalg.norm(g_leg[1][0:2, 2]-FL[0:2, 2])
    elif flag0 == 2: slip += np.linalg.norm(g_leg[2][0:2, 2]-HR[0:2, 2])
    elif flag0 == 3: slip += np.linalg.norm(g_leg[3][0:2, 2]-HL[0:2, 2])
    return slip


def func(xi_pre, l_a, l_b, l_c, l_d, l_e, alpha, beta, FR, FL, HL, HR, flag):
    target_func = lambda x: rms(np.sum(func1(x, l_c, beta, FR, FL, flag), 
                        axis = 1)) + rms(np.sum(func2(x, l_a, alpha, l_b, l_c,
                                                beta, HL, HR, flag), axis = 1))
    bnds = ((-np.inf, np.inf), (-np.inf, np.inf), (-np.inf, np.inf), 
            (l_d, l_d+l_e/4))
    # print(flag)
    # print("Function1: ", func1(xi_pre[0:4], l_c, beta, FR, FL, flag))
    # print("Function2: ", func2(xi_pre[0:4], l_a, alpha, l_b, l_c, 
    #                             beta, HL, HR, flag))
    # print("Initial value of target_func(1.2755e+01)", target_func(xi_pre[0:4]))
    res = minimize(target_func, xi_pre[0:4], bounds = bnds)
    res1 = fmin_slsqp(target_func, xi_pre[0:4], bounds = bnds, acc=1e-08)

    print("opt method1: ", res.x, "function value: ",)
    print("opt method2: ", res1)
    return res1

def func1(x, l_c, beta, FR, FL, flag):
    temp = G(x[0], x[1], x[2])
    if flag == 1:
        ans = matrix_multiply(temp, R(math.pi/2), F(l_c), R(-math.pi/2), 
                            R(beta[0]), F(x[3])) - FR
    else:
        ans = matrix_multiply(temp, R(-math.pi/2), F(l_c), R(math.pi/2), 
                            R(beta[1]), F(x[3])) - FL
    ans = ans*np.array([[0, 0, 1], [0, 0, 1],[0, 0, 0]])
    # print("func1: ", ans)
    return ans

def func2(x, l_a, alpha, l_b, l_c, beta, HL, HR, flag):
    temp = G(x[0], x[1], x[2])
    if flag == 1:
        ans = matrix_multiply(temp, F(l_a), R(alpha[0]), F(l_b), R(alpha[1]), 
                            F(l_a), R(-math.pi/2), F(l_c), R(math.pi/2), 
                            R(beta[3]), F(x[3])) - HL
    else:
        ans = matrix_multiply(temp, F(l_a), R(alpha[0]), F(l_b), R(alpha[1]), 
                            F(l_a), R(math.pi/2), F(l_c), R(-math.pi/2), 
                            R(beta[2]), F(x[3])) - HR
    ans = ans*np.array([[0, 0, 1], [0, 0, 1],[0, 0, 0]])
    # print("func2: ", ans)
    return ans