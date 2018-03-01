function [ gh, g_leg, xi, slip]=get_config(g_leg,beta,alpha,xi,leg_act,colorSpace)
load('robot_para.mat')

xi_pre=xi; %previous xi, also initialization
R=@(a)[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
F=@(L) [eye(2),[L;0];0 0 1];
FR=g_leg{1};FL=g_leg{2};HR=g_leg{3};HL=g_leg{4};
g=@(a,x,y) [cos(a) -sin(a) x;sin(a) cos(a) y;0 0 1];
air_leg_length=l_d+l_e/4;
slip=0;

% !!!!!!!!!!!!!!!!!!!!1
%load('g_leg_initial.mat');

myfun1_1=@(x) (g(x(1),x(2),x(3))*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(2))*F(x(4))-FL).*[0 0 1;0 0 1;0 0 0];
myfun1_2=@(x) (g(x(1),x(2),x(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(3))*F(x(4))-HR).*[0 0 1;0 0 1;0 0 0];
myfun2_1=@(x) (g(x(1),x(2),x(3))*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(2))*F(x(4))-FL).*[0 0 1;0 0 1;0 0 0];
myfun2_2=@(x) (g(x(1),x(2),x(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(3))*F(x(4))-HR).*[0 0 1;0 0 1;0 0 0];


if sum(leg_act) == 4
   myfun1 = @(x)abs(myfun1_1(x)) + abs(myfun1_2(x));
   myfun2 = @(x)abs(myfun2_1(x)) + abs(myfun2_2(x));

elseif leg_act(1)==1 && leg_act(4)==1
    myfun1=@(x) (g(x(1),x(2),x(3))*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(1))*F(x(4))-FR).*[0 0 1;0 0 1;0 0 0];
    myfun2=@(x) (g(x(1),x(2),x(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(4))*F(x(4))-HL).*[0 0 1;0 0 1;0 0 0];
elseif leg_act(2)==1 && leg_act(3)==1
    myfun1=@(x) (g(x(1),x(2),x(3))*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(2))*F(x(4))-FL).*[0 0 1;0 0 1;0 0 0];
    myfun2=@(x) (g(x(1),x(2),x(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(3))*F(x(4))-HR).*[0 0 1;0 0 1;0 0 0];
end
myfun=@(x) rms(sum(myfun1(x),2))+rms(sum(myfun2(x),2));
xi = fmincon(myfun,xi_pre(1:4),[],[],[],[],[-inf -inf -inf l_d],[inf inf inf l_d+l_e/4]);

%a = fmincon(myfun,x_test,[],[],[],[],[-inf -inf -inf l_d],[inf inf inf l_d+l_e/4]);

if sum(leg_act) == 4
    xi=[xi(1:4) xi(4) xi(4) xi(4)];
    slip=slip+norm(g_leg{2}(1:2,3)-FL(1:2,3)); % Dummy
    
elseif leg_act(1)==1 && leg_act(4)==1
    g_leg{1}=g(xi(1),xi(2),xi(3))*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(1))*F(xi(4));
    g_leg{2}=g(xi(1),xi(2),xi(3))*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(2))*F(air_leg_length);
    g_leg{3}=g(xi(1),xi(2),xi(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(3))*F(air_leg_length);
    g_leg{4}=g(xi(1),xi(2),xi(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(4))*F(xi(4));
    xi=[xi(1:4) air_leg_length air_leg_length xi(4)];
    if leg_act(2)==1
        slip=slip+norm(g_leg{2}(1:2,3)-FL(1:2,3));
    end
    if leg_act(3)==1
        slip=slip+norm(g_leg{3}(1:2,3)-HR(1:2,3));
    end
elseif leg_act(2)==1 && leg_act(3)==1
    g_leg{1}=g(xi(1),xi(2),xi(3))*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(1))*F(air_leg_length);
    g_leg{2}=g(xi(1),xi(2),xi(3))*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(2))*F(xi(4));
    g_leg{3}=g(xi(1),xi(2),xi(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(3))*F(xi(4));
    g_leg{4}=g(xi(1),xi(2),xi(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(4))*F(air_leg_length);
    xi=[xi(1:3) air_leg_length xi(4) xi(4) air_leg_length];
    if leg_act(1)==1
        slip=slip+norm(g_leg{1}(1:2,3)-FR(1:2,3));
    end
    if leg_act(4)==1
        slip=slip+norm(g_leg{4}(1:2,3)-HL(1:2,3));
    end
end


leg_act=[0 0 0 0 leg_act];

er=0.1;n=4;
gh=cell(1,5);
gh{1}=g(xi(1),xi(2),xi(3));
gh{2}=gh{1}*F(l_a)*R(alpha(1))*F(l_b/4);
gh{3}=gh{2}*F(l_b/2);
gh{4}=gh{3}*F(l_b/4)*R(alpha(2))*F(l_a);

%CoM=getCoM(gh,g_leg,leg_act,alpha(3));



end
