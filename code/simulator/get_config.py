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



def get_config(g_leg, betaFunc(params, t), xi, alphaFunc(params, t), blackred):
    
    sio.loadmat("robot_para.mat")
    # previous xi, also initialization
    xi_pre = xi

    

    return gh, g_leg, xi, slip    

function [h, gh, g_leg, xi, slip]=get_config(g_leg,beta,alpha,xi,leg_act,colorSpace)
load('robot_para.mat')

xi_pre=xi; %previous xi, also initialization
R=@(a)[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
F=@(L) [eye(2),[L;0];0 0 1];
FR=g_leg{1};FL=g_leg{2};HR=g_leg{3};HL=g_leg{4};
g=@(a,x,y) [cos(a) -sin(a) x;sin(a) cos(a) y;0 0 1];
air_leg_length=l_d+l_e/4;
slip=0;

if leg_act(1)==1 && leg_act(4)==1
    myfun1=@(x) (g(x(1),x(2),x(3))*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(1))*F(x(4))-FR).*[0 0 1;0 0 1;0 0 0];
    myfun2=@(x) (g(x(1),x(2),x(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(4))*F(x(4))-HL).*[0 0 1;0 0 1;0 0 0];
elseif leg_act(2)==1 && leg_act(3)==1
    myfun1=@(x) (g(x(1),x(2),x(3))*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(2))*F(x(4))-FL).*[0 0 1;0 0 1;0 0 0];
    myfun2=@(x) (g(x(1),x(2),x(3))*F(l_a)*R(alpha(1))*F(l_b)*R(alpha(2))*F(l_a)*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(3))*F(x(4))-HR).*[0 0 1;0 0 1;0 0 0];
end
myfun=@(x) rms(sum(myfun1(x),2))+rms(sum(myfun2(x),2));
xi = fmincon(myfun,xi_pre(1:4),[],[],[],[],[-inf -inf -inf l_d],[inf inf inf l_d+l_e/4]);

if leg_act(1)==1 && leg_act(4)==1
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

h=cell(1,13);

for i =[ 1 4 ]
    color = colorSpace(90+round(10*leg_act(i)),:);
    h{i} = drawActiveFrameEllipse(gh{i},l_a,er,color);
end
for i = [ 2 3]
    color = colorSpace(90+round(10*leg_act(i)),:);
    h{i} = drawActiveFrameEllipse(gh{i},l_b/4,er,color);
end
color = colorSpace(90+round(10*leg_act(i)),:);

h{13} = drawActiveFrameEllipse(gh{4}*F(l_b/4)*R(alpha(3))*F(l_a*4),l_a,er,color);

er=1;
color = colorSpace(90+round(54*leg_act(n+1)),:);
h{n+1} = drawActiveFrameEllipse(g_leg{1},l_a/5,er,color);
color = colorSpace(90+round(54*leg_act(n+2)),:);
h{n+2} = drawActiveFrameEllipse(g_leg{2},l_a/5,er,color);
color = colorSpace(90+round(54*leg_act(n+3)),:);
h{n+3} = drawActiveFrameEllipse(g_leg{3},l_a/5,er,color);
color = colorSpace(90+round(54*leg_act(n+4)),:);
h{n+4} = drawActiveFrameEllipse(g_leg{4},l_a/5,er,color);
sqs=[1 2;1 3;1 4;2 3;2 4;3 4];
line_ind=1;
for sqs_ind=1:6
    ind_1=sqs(sqs_ind,1)+4;
    ind_2=sqs(sqs_ind,2)+4;
    if leg_act(ind_1)*leg_act(ind_2)==1
        x1=g_leg{ind_1-4}(1:2,3);
        x2=g_leg{ind_2-4}(1:2,3);
        xx=[x1(1),x2(1)];
        yy=[x1(2),x2(2)];
        h{n+4+line_ind}=line(xx,yy);
        line_ind=line_ind+1;
    end
end
CoM=getCoM(gh,g_leg,leg_act,alpha(3));
beta=linspace(0,2*pi,101);
e_x=l_a/5*cos(beta);
e_y=l_a/5*sin(beta);
e_p=[e_x;e_y];

h{12} = patch(CoM(1)+e_p(1,:),CoM(2)+e_p(2,:),[0 0 1],'linewidth',1);

end

    return g_leg, xi, slip