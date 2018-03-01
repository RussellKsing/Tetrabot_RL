function CoM=getCoM(gh,g_leg,leg_act,alpha)
load('robot_para.mat')
F=@(L) [eye(2),[L;0];0 0 1];
R=@(a)[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];

for ind=1:4
    if leg_act(ind+4)==0
        g_leg{ind}=g_leg{ind}*F(l_e*3/4);
    end
end

CoM=[0;0];
% for ind=1:4
%     CoM=CoM+1*gh{ind}(1:2,3);
%     CoM=CoM+1*g_leg{ind}(1:2,3);
% end
gh5=gh{4}*F(l_b/4)*R(alpha)*F(l_a*4);
CoM=CoM+3*gh5(1:2,3);
CoM=CoM+3*gh{1}(1:2,3);
CoM=CoM+1*gh{2}(1:2,3);
CoM=CoM/7;