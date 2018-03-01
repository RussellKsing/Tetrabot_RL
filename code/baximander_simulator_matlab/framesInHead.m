function [ g, gFL, gFR, gHL, gHR ] = framesInHead( alpha, beta, hind)
n=size(alpha,1)+1;
g=cell(1,n);
load('robot_para.mat')
g{1}=eye(3);
F=@(L) [eye(2),[L;0];0 0 1];%move forward
R=@(a)[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];%rotation frame
g{2}=g{1}*F(l_a)*R(alpha(1))*F(l_b/4);
g{3}=g{2}*F(l_b/2);
g{4}=g{3}*F(l_b/4)*R(alpha(2))*F(l_a);

gFL = g{1}*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(1))*[eye(2),[l_d+l_e/4;0];0 0 1];
gFR = g{1}*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(2))*[eye(2),[l_d+l_e/4;0];0 0 1];
gHL = g{hind-1}*R(pi/2)*F(l_c)*R(-pi/2)*R(beta(3))*[eye(2),[l_d+l_e/4;0];0 0 1];
gHR = g{hind-1}*R(-pi/2)*F(l_c)*R(pi/2)*R(beta(4))*[eye(2),[l_d+l_e/4;0];0 0 1];

end

