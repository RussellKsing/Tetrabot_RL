# main page for salamander
from bot import *
import scipy.io
# global variables
# dutyf, interPhi, meterPhi, outerPhi


def main(params):
    blackred = sio.loadmat("mat_files/BlackWhiteRedColormap.mat")

    Aleg = math.pi/3
   
    # current time
    t_ini = t = 0
    # Total number of time steps
    T = 100
    # time period
    dt = 2*math.pi/T
    f_num = 0
    alpha_pre_col = leg_col = beta_col = []
    
    # this could generate an error: l_d, l_e not initialized ??
    xi = np.array([-math.pi/2, 0, 0, l_d+l_e/2, l_d+l_e/2, 
                l_d+l_e/2, l_d+l_e/2])
    h_ini = [0, 0]

    bodyContact = [1.2, 1.2, 1.2, 0]
    leg_act = [1,1,1,1]
    init_act = np.array(bodyContact + leg_act)

    g_h = np.array([[math.cos(q_h[2]), -math.sin(q_h[2]), q_h[0],],
                    [math.sin(q_h[2]),  math.cos(q_h[2]), q_h[1]],
                    [               0,                 0,      1]])

    # [ h ] = drawActiveFrameSnake( [0 1 h_ini(1);-1 0 h_ini(2);0 0 1], alpha(t), beta(t), 5, [0 0 0 0 1 1 1 1], blackred );

    h = drawActiveFrameSnake(g_h, alphaFunc(params, t), betaFunc(params, t),
                                 hind, L, Lleg, init_act, blackred)


    data = scipy.io.loadmat('g_leg_initial.mat')
    g_leg = np.array(data['g_leg'])

    FR_length=np.zeros(1,100); FR_ang=np.zeros(1,100);
    FL_length=np.zeros(1,100); FL_ang=np.zeros(1,100);
    HR_length=np.zeros(1,100); HR_ang=np.zeros(1,100);
    HL_length=np.zeros(1,100); HL_ang=np.zeros(1,100);
    body_length=np.zeros(1,100);body_lateral=np.zeros(1,100);
    slip_col=0;


    while(t < t_ini + 2*math.pi):
        g_leg, xi, slip = get_config(g_leg, betaFunc(params, t), xi, alphaFunc(params, t), leg_act, blackred)
        slip_col = slip_col + slip
        act_cur = activationFunc(bodyContact, params, t)
        act_next= activationFunc(bodyContact, params, t+dt)

        # drawnow

        f_num = f_num+1;
        t = t + dt

        FR_length(f_num) = xi[3] * f_leg_act(dutyf, t)
        FR_ang[f_num] = f_leg(dutyf, Aleg, t)

        FL_length(f_num) = xi[4] * f_leg_act(dutyf, t)
        FL_ang[f_num] = f_leg(dutyf, Aleg, t)

        HR_length(f_num) = xi[5] * f_leg_act(dutyf, t)
        HR_ang[f_num] = f_leg(dutyf, Aleg, t)

        HL_length(f_num) = xi[6] * f_leg_act(dutyf, t)
        HL_ang[f_num] = f_leg(dutyf, Aleg, t)

        # ??
        body_length[f_num] = l_b
        body_lateral[f_num] = params.a2 * math.sin(t + params.b2)



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
# beta = @(c)[c_FR*F_Leg(Aleg,c,dutyf), -c_FL*F_Leg(Aleg,c+interPhi,dutyf), c_HR*F_Leg(Aleg,c+meterPhi,dutyf), -c_HL*F_Leg(Aleg,c+outerPhi,dutyf)] + beta0;
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
    return np.array([fct(params, t, 0), fct(params, t, 1), 0.])

# this function computes the derivative of the alpha variable
def d_alphaFunc(params, t):
    return np.array([(fct(params, t + 0.01, 0) - fct(params, t, 0))/0.01, (fct(params, t + 0.01, 1) - fct(params, t, 1))/0.01, 0.])

def fct(params, t, i):
    '''
    IF 9params, Fct = @(a1,a2,b1,b2,b3) (a1*sin(t/2+b1)+a2*sin(t+b2)+a3*sin(t*2+b3));
    :param params:
        all constants
        params.a1, params.a2, params.b1, params.b2, params.b3
    '''
    b1 = 0.
    return params.a1[i] * math.sin(t / 2. + params.b1[i]) + params.a2[i] * math.sin(t + params.b2[i])

##### Draw Active Snake ######
# this function returns the new position of the salamander given the
# previous motion and the equation of motion given by alpha, beta.
def drawActiveFrameSnake(g, alpha, beta, hind, L, Lleg, activation, colorSpace):
    n = alpha.shape[0] + 1
    [k, gFR, gFL, gHR, gHL] = framesInHead(alpha, beta, hind, L, Lleg)

    gFR = np.dot(g, gFR,)
    gFL = np.dot(g, gFL)
    gHL = np.dot(g, gHL)
    gHR = np.dot(g, gHR)

    for i in xrange(0, n):
        k[i] = np.dot(g, k[i])

    massCenter, linksCenter = computeCOM(k, activation[:4])

    endPoints = computeEND(k, L, activation[:4])
    
    posFL = (gFL[0, 2], gFL[1, 2])
    posFR = (gFR[0, 2], gFR[1, 2])
    posHL = (gHL[0, 2], gHL[1, 2])
    posHR = (gHR[0, 2], gHR[1, 2])
    posMass = (massCenter[0, 0], massCenter[0, 1])

    # activation = [bodyx4, front right, front left, back right, back left]
    curPos = (posFR, posFL, posHL, posHR)
    curAct = (activation[4], activation[5], activation[7], activation[6])

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

    stablePenaltyScore = checkStability(legsOnGround_x, legsOnGround_y, posMass)

    # pprint (stablePenaltyScore)
    # sys.exit()
    # For plotting a triangle for the legs on the ground
    if (len(legsOnGround_x) != 0):
        legsOnGround_x += [legsOnGround_x[0]]
        legsOnGround_y += [legsOnGround_y[0]]

    newPoints = (legsOnGround_x, legsOnGround_y, legsInAir_x, legsInAir_y, posMass[0], posMass[1])

    # Return
    return (massCenter, newPoints, stablePenaltyScore, linksCenter, endPoints)

def computeCOM(g, activation):
    com = np.zeros([1,2], dtype=float)
    numlinks = np.count_nonzero(activation)
    comlist = []
    # print (numlinks)
    for i in xrange(0, numlinks):
        # print (g[i])
        # print ('inside', g[i, 0:2, 2])
        com = com + g[i, 0:2, 2]
        comlist.append(g[i, 0:2, 2])
        # print ('com', com)
    xcom = [item[0] for item in comlist]
    ycom = [item[1] for item in comlist]
    return com/numlinks, [xcom, ycom]

def computeEND(k, L, activation):
    numlinks = np.count_nonzero(activation)
    F = np.array([[1,0,L],[0,1,0],[0,0,1]])
    comlist = []
    # end_matrix = []
    end_matrix = np.zeros((4,3,3))
    for i in xrange(0, numlinks):
        end_matrix[i] = np.dot(k[i], np.linalg.inv(F))

    end_matrix[numlinks] = np.dot(k[i], F)

    for i in xrange(0, numlinks+1):
        comlist.append(end_matrix[i, 0:2, 2])
    
    x_endlist = [item[0] for item in comlist]
    y_endlist = [item[1] for item in comlist]

    return [x_endlist, y_endlist]

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

# forward matrix
def F(L):
    ans = np.array([[1, 0, L], [0, 1, 0], [0, 0, 1]])
    return ans

# rotation matrix
# a = angle to rotate by
def R(a):
    ans = np.array([[math.cos(a), -math.sin(a), 0],
           [math.sin(a), math.cos(a), 0],
           [0, 0, 1]])
    return ans 

def G(a, x, y):
    ans = np.array([[math.cos(a), -math.sin(a), x],
         [math.sin(a), math.cos(a), y], [0, 0, 1]])
    return ans

def get_config(g_leg, beta, xi, alpha, leg_act, blackred):
    
    sio.loadmat("robot_para.mat")
    # previous xi, also initialization
    xi_pre = xi
    
    FR = g_leg[0]
    FL = g_leg[1]
    HR = g_leg[2]
    HL = g_leg[3]

    air_leg_length = (l_d + l_e)/4
    slip = 0

    xi = fmincon(myfun,xi_pre(1:4),[],[],[],[],[-inf -inf -inf l_d],[inf inf inf l_d+l_e/4]);

    if leg_act[0]==1 and leg_act[4]==1:
        g_leg[0] = np.dot(np.dot(np.dot(np.dot(np.dot(G(xi[0], xi[1], xi[2]), R(math.pi/2)), F(l_c)), R(-math.pi/2)), R(beta[0])), F(xi[3]))
        g_leg[1] = np.dot(np.dot(np.dot(np.dot(np.dot(G(xi[0], xi[1], xi[2]), R(-math.pi/2)), F(l_c)), R(math.pi/2)), R(beta[1])), F(air_leg_length))
        g_leg[2] = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(G(xi[0], xi[1], xi[2]), F(l_a)), R(alpha[0])), F(l_b)), R(alpha[1])), F(l_a)), R(math.pi/2)), F(l_c)), R(-math.pi/2)), R(beta[2])), F(air_leg_length))
        g_leg[3] = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(G(xi[0], xi[1], xi[2]), F(l_a)), R(alpha[0])), F(l_b)), R(alpha[1])), F(l_a)), R(-math.pi/2)), F(l_c)), R(math.pi/2)), R(beta[3])), F(xi[3]))
       
        xi = [xi[:3], air_leg_length, air_leg_length, xi[3]]

        if leg_act[1] == 1:
            slip = slip + norm(g_leg[1][:2, 2] - FL[:2, 3])
        if leg_act[2] == 1:
            slip = slip + norm(g_leg[1][:2, 2] - HR[:2, 3])
    
    elif leg_act[1]==1 and leg_act[2]==1:
        g_leg[0] = np.dot(np.dot(np.dot(np.dot(np.dot(G(xi[0], xi[1], xi[2]), R(math.pi/2)), F(l_c)), R(-math.pi/2)), R(beta[0])), F(air_leg_length))
        g_leg[1] = np.dot(np.dot(np.dot(np.dot(np.dot(G(xi[0], xi[1], xi[2]), R(-math.pi/2)), F(l_c)), R(math.pi/2)), R(beta[1])), F(xi[3]))
        g_leg[2] = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(G(xi[0], xi[1], xi[2]), F(l_a)), R(alpha[0])), F(l_b)), R(alpha[1])), F(l_a)), R(math.pi/2)), F(l_c)), R(-math.pi/2)), R(beta[2])), F(xi[3]))
        g_leg[3] = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(G(xi[0], xi[1], xi[2]), F(l_a)), R(alpha[0])), F(l_b)), R(alpha[1])), F(l_a)), R(-math.pi/2)), F(l_c)), R(math.pi/2)), R(beta[3])), F(air_leg_length))
       
        xi = [xi[:3], air_leg_length, xi[3], xi[3], air_leg_length]

        if leg_act[0] == 1:
            slip = slip + norm(g_leg[1][:2, 2] - FR[:2, 3])
        if leg_act[3] == 1:
            slip = slip + norm(g_leg[1][:2, 2] - HR[:2, 3])

        leg_act = [0, 0, 0, 0] + leg_act
        er=0.1; n=4;
        gf = cell[0,5]

        gh[0] = G(xi[0],xi[1],xi[2]);
        gh[1] = np.dot(np.dot(np.dot(gh[0], F(l_a)), R(alpha[0])), F(l_b/4))
        gh[2] = np.dot(gh[1], F(l_b/2))
        gh[3] = np.dot(np.dot(np.dot(np.dot(gh[2], F(l_b/4))), R(alpha[1])), F(l_a))

        # ??
        h = cell[0, 12]

    return gh, g_leg, xi, slip

if __name__ == '__main__':
    class Struct(object): pass
    params = Struct()
    params.numParams = 11 # if != 9, then fct return 0
    params.Aleg = math.pi / 3
    params.slope = 0

    params.dutyf = dutyFactor
    params.interPhi = interPhi
    params.meterPhi = meterPhi
    params.outerPhi = outerPhi

    params.a2 = [a2_1, a2_2]
    params.b2 = [b2_1, b2_2]

    params.b1 = [0., 0.]
    params.a1 = [0., 0.]
    params.a3 = [0., 0.]

    # turnning scaler for beta (leg angle)
    params.cFL = c1
    params.cHL = c2
    params.cHR = c3
    params.cFR = c4

    ######## CHANGE IT TO FALSE TO TURN OFF PLOTTING ########
    params.ifPlot = True


    # Start to run the simulator with different duty_ind
    (angleChange, displacement, unstablePenalty) = main(params)
    if (printResult):
        print "Angle Change: " + str(angleChange) + " pi"
        print "Displacement: " + str(displacement) + "cm"
        print "Unstable Penalty: " + str(unstablePenalty) + "\n"

    return (angleChange, displacement, unstablePenalty)
