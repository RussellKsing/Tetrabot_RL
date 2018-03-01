%% Initialization
Aleg=pi/3;


leg_act =@(t) [F_leg_act(t,dutyf), F_leg_act(t+interPhi,dutyf), F_leg_act(t+meterPhi,dutyf), F_leg_act(t+outerPhi,dutyf)];
beta0=[1 -1 1 -1]*pi/2;
beta = @(c)[c_FR*F_Leg(Aleg,c,dutyf), -c_FL*F_Leg(Aleg,c+interPhi,dutyf), c_HR*F_Leg(Aleg,c+meterPhi,dutyf), -c_HL*F_Leg(Aleg,c+outerPhi,dutyf)] + beta0;
%alpha=@(t) [0;0];
%current time
t_ini=0;
T = 100;%number of time steps
dt = 2*pi/T;

t = t_ini; %!!!!!!!!!!!!!!!!!

xi=[-pi/2 0 0 l_d+l_e/2 l_d+l_e/2 l_d+l_e/2 l_d+l_e/2];
h_ini=[0;0]; % head initial position


load('g_leg_initial.mat')    
%% Simulation
%for ind=1:5
blackred = 1; % Dummy value
i = 0;

g_leg_init = g_leg;

while t<t_ini+2*pi  
    i = i + 1;
    % g_leg = g_leg_init;
    [gh, g_leg,xi, slip]=get_config(g_leg,beta(t),alpha(t),xi,leg_act(t),blackred);

    t=t+dt;
end
%load('robot_para.mat')
gh{1}(1:2,3)-h_ini
norm(gh{1}(1:2,3)-h_ini)
