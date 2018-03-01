function [ h,Rs,g,g_leg] = drawsalamander( leg_act, g_leg, beta, L, Lleg, colorSpace,alpha_pre )

option=optimset('Display','notify');
g=cell(1,4);
R=@(a)[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
F=@(L) [eye(2),[L;0];0 0 1];
Lsqaure_1=3.1;
Lsqaure_2=3.6;
if leg_act(1)==1 && leg_act(2)==1
    gFL=g_leg{2};gFR=g_leg{1};
    g{1}=gFL/F(Lleg)/R(beta(2))/R(pi/2)/F(Lsqaure_1)/R(-pi/2);
    if leg_act(3)==1
        myfun=@(Rs) (g{1}*F(Lsqaure_2)*R(Rs(1))*F(2*L)*R(alpha_pre(3))*F(2*L)*R(Rs(2))*F(Lsqaure_2)*R(pi/2)*F(Lsqaure_1)*R(-pi/2)*R(beta(3))*F(Lleg));
        myfun1=@(Rs) sum(rms((myfun(Rs)-g_leg{3}).*[0 0 1;0 0 1;0 0 0]));
        Rs=fminsearch(@(Rs) myfun1(Rs),[alpha_pre],option);
    else
        myfun=@(Rs) (g{1}*F(Lsqaure_2)*R(Rs(1))*F(2*L)*R(alpha_pre(3))*F(2*L)*R(Rs(2))*F(Lsqaure_2)*R(-pi/2)*F(Lsqaure_1)*R(pi/2)*R(beta(4))*F(Lleg));
        myfun1=@(Rs) sum(rms((myfun(Rs)-g_leg{4}).*[0 0 1;0 0 1;0 0 0]));
        Rs=fminsearch(@(Rs) myfun1(Rs),[alpha_pre],option);        
    end
    g{2}=g{1}*F(Lsqaure_2)*R(Rs(1))*F(L);
    g{3}=g{2}*F(L)*R(alpha_pre(3))*F(L);
    g{4}=g{3}*F(L)*R(Rs(2))*F(Lsqaure_2);
    gHR=g{4}*R(pi/2)*F(Lsqaure_1)*R(-pi/2)*R(beta(3))*F(Lleg);
    gHL=g{4}*R(-pi/2)*F(Lsqaure_1)*R(pi/2)*R(beta(4))*F(Lleg);
else
    gHL=g_leg{4};gHR=g_leg{3};
    g{4}=gHL/F(Lleg)/R(beta(4))/R(pi/2)/F(Lsqaure_1)/R(-pi/2);
    if leg_act(1)==1
        myfun=@(Rs) (g{4}/F(Lsqaure_2)/R(Rs(2))/F(2*L)/R(alpha_pre(3))/F(2*L)/R(Rs(1))/F(Lsqaure_2)*R(pi/2)*F(Lsqaure_1)*R(-pi/2)*R(beta(1))*F(Lleg));
        myfun1=@(Rs) sum(rms((myfun(Rs)-g_leg{1}).*[0 0 1;0 0 1;0 0 0]));
        Rs=fminsearch(@(Rs) myfun1(Rs),[alpha_pre],option);
    else
        myfun=@(Rs) (g{4}/F(Lsqaure_2)/R(Rs(2))/F(2*L)/R(alpha_pre(3))/F(2*L)/R(Rs(1))/F(Lsqaure_2)*R(-pi/2)*F(Lsqaure_1)*R(pi/2)*R(beta(2))*F(Lleg));
        myfun1=@(Rs) sum(rms((myfun(Rs)-g_leg{2}).*[0 0 1;0 0 1;0 0 0]));
        Rs=fminsearch(@(Rs) myfun1(Rs),[alpha_pre],option);     
    end
    g{3}=g{4}/F(Lsqaure_2)/R(Rs(2))/F(L);
    g{2}=g{3}/F(L)/R(alpha_pre(3))/F(L);
    g{1}=g{2}/F(L)/R(Rs(1))/F(Lsqaure_2);
    gFR=g{1}*R(pi/2)*F(Lsqaure_1)*R(-pi/2)*R(beta(1))*F(Lleg);
    gFL=g{1}*R(-pi/2)*F(Lsqaure_1)*R(pi/2)*R(beta(2))*F(Lleg);
end
leg_act=[0 0 0 0 leg_act];
er=0.1;n=4;
for i =[ 1 4 ]
    color = colorSpace(90+round(10*leg_act(i)),:);
    h{i} = drawActiveFrameEllipse(g{i},2*L,er,color);
end
for i = [ 2 3]
    color = colorSpace(90+round(10*leg_act(i)),:);
    h{i} = drawActiveFrameEllipse(g{i},1.5*L,er,color);
end
er=1;
color = colorSpace(90+round(54*leg_act(n+1)),:);
h{n+1} = drawActiveFrameEllipse(gFR,L/5,er,color);
color = colorSpace(90+round(54*leg_act(n+2)),:);
h{n+2} = drawActiveFrameEllipse(gFL,L/5,er,color);
color = colorSpace(90+round(54*leg_act(n+3)),:);
h{n+3} = drawActiveFrameEllipse(gHR,L/5,er,color);
color = colorSpace(90+round(54*leg_act(n+4)),:);
h{n+4} = drawActiveFrameEllipse(gHL,L/5,er,color);

g_leg{1}=gFR;
g_leg{2}=gFL;
g_leg{3}=gHR;
g_leg{4}=gHL;
Rs=Rs(1:2);
