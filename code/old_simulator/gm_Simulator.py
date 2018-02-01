import numpy as np
import math
import scipy.io as sio
from scipy.linalg import expm, norm
from computeBodyVelocity import *
import matplotlib.pyplot as plt
import sys
from pprint import pprint

def gm_Simulator(params):
    blackred = sio.loadmat("mat_files/BlackWhiteRedColormap.mat")
    N = 4  # number of modules
    undulationNumber = 1
    phaseOffset = math.pi

    L = 5.5 # L is body length
    Lleg = 7 # Lleg is leg length
    hind = 3

    t = np.array([0.]) # current time
    T = 300 # number of time steps
    dt = 2.0 * math.pi/T


    # beta = betaFunc(params, params.c)
    # d_beta = d_betaFunc(params, params.c)
    # alpha = alphaFunc(params, t)
    # d_alpha = d_alphaFunc(params, t)

    joint_index = np.array(range(1, N)) # joint index

    # what is bodycontact?
    bodyContact = [1.2, 1.2, 0, 0]

    # what is init_act?
    init_act = np.array(bodyContact + [1, 1, 1, 1])

    # Initial head position, [x,y,theta]
    q_h = np.array([0, 0, -math.pi/2.])


    g_h = np.array([[math.cos(q_h[2]), -math.sin(q_h[2]), q_h[0],],
                    [math.sin(q_h[2]),  math.cos(q_h[2]), q_h[1]],
                    [               0,                 0,      1]])

    # print (alphaFunc(params, t))
    # print (betaFunc(params, t))
    # sys.exit()
    # comInit - initial center of mass
    # newPoints - legsOnGround_x, legsOnGround_y, legsInAir_x, legsInAir_y, com_x, com_y
    (comInit, newPoints, stablePenaltyScore) = drawActiveFrameSnake(g_h, alphaFunc(params, t), betaFunc(params, t),
                                                                    hind, L, Lleg, init_act, blackred)
    # pprint (comInit)
    # pprint (newPoints)
    # pprint (stablePenaltyScore)
    # sys.exit()

    # friction profile
    K = np.diag([1,2,1]) 

    F = np.array([-4.5 * math.sin(params.slope * math.pi/180.), 0, 0])
    # pprint(F)
    # sys.exit()
    ## Init plot
    if (params.ifPlot):
        fig, ax = plt.subplots()
        legsOnGround, = ax.plot(newPoints[0], newPoints[1], marker='o', color='b', linestyle='-')
        legsInTheAir, = ax.plot(newPoints[2], newPoints[3], marker='o', color='y', linestyle='None')
        bodyMassPos, = ax.plot(newPoints[4], newPoints[5], marker='o', color='r', linestyle='None')
        ax.set_xlim(-40, 40)
        ax.set_ylim(-40, 40)
        plt.pause(6)


    unstablePenalty = stablePenaltyScore
    itr = 1
    while (t[0] < math.pi * 2 + 0.01):
        itr += 1
        #print t
        activation1 = activationFunc(bodyContact, params, t)
        # pprint (activation1)
        # sys.exit()

        # xi has three components
        xi = computeBodyVelocity(alphaFunc(params, t), d_alphaFunc(params, t),
                                 betaFunc(params, t), d_betaFunc(params, t),
                                 hind, L, Lleg, activation1, K, F)

        # pprint (xi)
        # sys.exit()

        if (itr % 1 == 0):
            (curCom, newPoints, stablePenaltyScore) = drawActiveFrameSnake(g_h, alphaFunc(params, t), betaFunc(params, t),
                                                                            hind, L, Lleg, activation1, blackred)

            # Accumulate the penalty score
            #print stablePenaltyScore
            #print unstablePenalty
            unstablePenalty += stablePenaltyScore


        # Update plots
        if (params.ifPlot):
            legsOnGround.set_data(newPoints[0], newPoints[1])
            legsInTheAir.set_data(newPoints[2], newPoints[3])
            bodyMassPos.set_data(newPoints[4], newPoints[5])
            plt.pause(0.015)


        t += dt
        g_h = np.dot(g_h, expm(dt * xiHat(xi))) # move the head
        #print xi
        #print g_h


    itr += 1.
    (comEnd, newPoints, stablePenaltyScore) = drawActiveFrameSnake(g_h, alphaFunc(params, t), betaFunc(params, t),
                                                                   hind, L, Lleg, init_act, blackred)
    unstablePenalty += stablePenaltyScore


    # angle = 0.5 + 1 / math.pi * math.asin(-g_h(1,2))

    angleChange = 0.5 + 1./math.pi * np.arcsin(-g_h[0,1])
    displacement = norm(comEnd - comInit)

    unstableScaler = unstablePenalty / (itr)

    return (angleChange, displacement, unstableScaler)



def xiHat(xi):
    return np.array([[0, -xi[2], xi[0],],
                     [xi[2], 0, xi[1]],
                     [    0, 0,     0]])



##### Activation function ########
def activationFunc(bodyContact, params, t):
    return np.concatenate([np.array(bodyContact),
                           f_leg_act(params.dutyf, t),
                           f_leg_act(params.dutyf, t + params.interPhi),
                           f_leg_act(params.dutyf, t + params.meterPhi),
                           f_leg_act(params.dutyf, t + params.outerPhi)])

# function to compute leg activation - whether the leg is on the ground or not?
def f_leg_act(dutyf, time):
    leg_act_col = np.zeros(time.shape)
    # pprint (dutyf)
    # pprint (time)
    for t_ind in xrange(0, time.size):
        t = time[t_ind]
        while (t < 0 or t > 2 * math.pi):
            if (t < 0):
                t = t + 2 * math.pi
            elif (t > 2 * math.pi):
                t = t - 2 * math.pi
        if (t < (1 - dutyf) * math.pi):
            leg_act = 0
        elif (t < 2 * math.pi - (1 - dutyf) * math.pi):
            leg_act = 1
        else:
            leg_act = 0
        leg_act_col[t_ind] = leg_act
    return leg_act_col



##### Beta function ######
def betaFunc(params, c):
    beta0 = np.array([math.pi / 2, -math.pi / 2., math.pi / 2., -math.pi / 2.], dtype=np.float32)
    dutyf = params.dutyf
    Aleg = params.Aleg
    beta =  np.array([(f_leg(dutyf, Aleg, c))                     * params.cFR,
                      (-f_leg(dutyf, Aleg, c + params.interPhi))  * params.cFL,
                      (f_leg(dutyf, Aleg, c + params.meterPhi))   * params.cHR,
                      (-f_leg(dutyf, Aleg, c + params.outerPhi))  * params.cHL])
    # print (beta)
    # print (np.reshape(beta, [1,4])[0]) + beta0
    # sys.exit()
    return np.reshape(beta, [1,4])[0] + beta0

def d_betaFunc(params, c):
    # d_beta = @(c)(beta(c+0.01)-beta(c))/0.01;
    return (betaFunc(params, c + 0.01) - betaFunc(params, c))/0.01

# what does this function do?
def f_leg(dutyf, Aleg, time):
    leg_angle_col = np.zeros(time.shape)
    for t_ind in xrange(0, time.size):
        t = time[t_ind]
        # conversion of range - to between 0 to 2pi
        while (t < 0 or t > 2*math.pi):
            if (t < 0):
                t = t + 2 * math.pi
            elif (t > 2 * math.pi):
                t = t - 2 * math.pi

        if (t < (1- dutyf) * math.pi):
            leg_angle = Aleg * math.sin(t / 2. / (1 - dutyf))
        elif ( t < 2 * math.pi - (1-dutyf) * math.pi):
            leg_angle = Aleg * math.cos((t-(1-dutyf)*math.pi)*(math.pi/(2*math.pi-2*(1-dutyf)*math.pi)))
        else:
            leg_angle = -Aleg * math.cos((t-2*math.pi+(1-dutyf)*math.pi)/2./(1-dutyf))

        leg_angle_col[t_ind] = leg_angle
        ttt = leg_angle_col[t_ind]
    return leg_angle_col


##### Alpha function ######
def alphaFunc(params, t):
    # [bodyx4, front left, front right, back left, back right]
    return np.array([fct(params, t), 0.,0.])

# this function computes the derivative of the alpha variable
def d_alphaFunc(params, t):
    return np.array([(fct(params, t + 0.01) - fct(params, t))/0.01, 0., 0.])

def fct(params, t):
    '''
    IF 9params, Fct = @(a1,a2,b1,b2,b3) (a1*sin(t/2+b1)+a2*sin(t+b2)+a3*sin(t*2+b3));
    :param params:
        all constants
        params.a1, params.a2, params.b1, params.b2, params.b3
    '''
    b1 = 0.
    return params.a1 * math.sin(t / 2. + b1) + params.a2 * math.sin(t + params.b2)


'''
class Struct(object):pass
params = Struct()
params.numParams = 9
params.a1 = 0
params.a2 = math.pi/3.
params.b1 = 0.
params.b2 = math.pi
params.b3 = 0
print fct(params, 1)
'''

##### Draw Active Snake ######
# this function returns the new position of the salamander given the
# previous motion and the equation of motion given by alpha, beta.
def drawActiveFrameSnake(g, alpha, beta, hind, L, Lleg, activation, colorSpace):
    n = alpha.shape[0] + 1
    [k, gFR, gFL, gHR, gHL] = framesInHead(alpha, beta, hind, L, Lleg)
    # pprint (g)
    # pprint (k)
    # pprint (gFR)
    # pprint (gFL)
    # pprint (gHR)
    # pprint (gHL)
    # sys.exit()

    gFR = np.dot(g, gFR,)
    gFL = np.dot(g, gFL)
    gHL = np.dot(g, gHL)
    gHR = np.dot(g, gHR)

    # pprint (g)
    # pprint (k)
    # pprint (gFR)
    # pprint (gFL[0, 2])#, gFL[1, 2])
    # pprint (gHR)
    # pprint (gHL)
    # sys.exit()

    k[0] = g
    for i in xrange(1, n-2):
        k[i] = np.dot(g, k[i])
    k[n-1] = np.dot(g, k[n-1])
    massCenter = computeCOM(k)

    # print (massCenter[0])
    # sys.exit()
    posFL = (gFL[0, 2], gFL[1, 2])
    posFR = (gFR[0, 2], gFR[1, 2])
    posHL = (gHL[0, 2], gHL[1, 2])
    posHR = (gHR[0, 2], gHR[1, 2])
    posMass = (massCenter[0, 0], massCenter[0, 1])

    # activation = [bodyx4, front right, front left, back right, back left]
    curPos = (posFR, posFL, posHL, posHR)
    curAct = (activation[4], activation[5], activation[7], activation[6])

    # pprint (curPos)
    # pprint (curAct)
    # sys.exit()
    # Get the list of legs that on the ground / in the air
    # legsOnGround_x - the x coordinate of the legs on ground
    # legsInAir_x - previous x coordinate of the legs in the air?
    (legsOnGround_x, legsOnGround_y, legsInAir_x, legsInAir_y) = ([], [], [], [])
    for i in range(4):
        if (curAct[i] == 1):
            legsOnGround_x += [curPos[i][0]]
            legsOnGround_y += [curPos[i][1]]
        else:
            legsInAir_x += [curPos[i][0]]
            legsInAir_y += [curPos[i][1]]

    # pprint (legsOnGround_x)
    # pprint (legsOnGround_y)
    # pprint (legsInAir_x)
    # pprint (legsInAir_y)
    # sys.exit()
    # Check stability
    stablePenaltyScore = checkStability(legsOnGround_x, legsOnGround_y, posMass)

    # pprint (stablePenaltyScore)
    # sys.exit()
    # For plotting a triangle for the legs on the ground
    if (len(legsOnGround_x) != 0):
        legsOnGround_x += [legsOnGround_x[0]]
        legsOnGround_y += [legsOnGround_y[0]]

    newPoints = (legsOnGround_x, legsOnGround_y, legsInAir_x, legsInAir_y, posMass[0], posMass[1])

    # Return
    return (massCenter, newPoints, stablePenaltyScore)

def computeCOM(g):
    com = np.zeros([1,2], dtype=float)
    for i in xrange(0, 2):
        com = com + g[i, 0:2, 2]
    return com/2.

def checkStability(legsOnGround_x, legsOnGround_y, posMass):
    # Given a score
    # If only 0, 1, 2 legs on the ground, is unstable; if all legs on the ground, stable
    if (len(legsOnGround_x) == 0 or len(legsOnGround_x) == 1 or len(legsOnGround_x) == 2):
        return 1.

    elif (len(legsOnGround_x) == 4):
        return 0.

    # 4 legs on the ground
    (x1, y1) = (legsOnGround_x[0], legsOnGround_y[0])
    (x2, y2) = (legsOnGround_x[1], legsOnGround_y[1])
    (x3, y3) = (legsOnGround_x[2], legsOnGround_y[2])

    # Check stability
    (inBound1, distance1) = checkInBound(x1, y1, x2, y2, posMass[0], posMass[1])
    (inBound2, distance2) = checkInBound(x2, y2, x3, y3, posMass[0], posMass[1])
    (inBound3, distance3) = checkInBound(x3, y3, x1, y1, posMass[0], posMass[1])

    if (inBound1 < 0 and inBound2 < 0 and inBound3 < 0):
        return 0.
    else:
        return 1.

def checkInBound(x1, y1, x2, y2, checkX, checkY):
    '''
    If (x, y) within bound, the point is in the tra
    test1 = @(x, y) (x-x1)*(y2-y1) - (y - y1)*(x2-x1);
    test2 = @(x, y) (x-x2)*(y3-y2) - (y - y2)*(x3-x2);
    test3 = @(x, y) (x-x3)*(y1-y3) - (y - y3)*(x1-x3);
    '''
    result = (checkX - x1) * (y2 - y1) - (checkY - y1) * (x2 - x1)
    distance = distancePointToLine(x1, y1, x2, y2, checkX, checkY)
    #print 'result:'
    #print result
    return (result, distance)

def distancePointToLine(x1, y1, x2, y2, xCheck, yCheck):
    line = abs((xCheck - x1) * (y2 - y1) - (yCheck - y1) * (x2 - x1))
    return line/math.sqrt(math.pow(y2-y1, 2) + math.pow(x2-x1, 2))

def getPenaltyScore(distance, acceptableBound):
    return -pow(2, (3 * (distance - acceptableBound)))/100.


'''
class Struct(object): pass
params = Struct()
params.Am = 3 * math.pi / 12
params.Aleg = math.pi / 3
params.interPhi = 3.1416
params.meterPhi = 0.3142
params.outerPhi = -2.8274
params.dutyf = 0.9
bodyContact = [1.2, 1.2, 0, 0]
t = np.array([5.6549])
'''



