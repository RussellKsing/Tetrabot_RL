function [ h ] = drawActiveFrameSnake( g, alpha, beta, hind, activation, colorSpace )
load('robot_para.mat')
n = size(alpha,1)+1;
h = cell(1,13);
color = colorSpace(90+round(10*activation(1)),:);
er=0.1;
h{1} = drawActiveFrameEllipse(g, l_a,er, color); %return
[k, gFL, gFR, gHL, gHR] = framesInHead(alpha, beta, hind);
k{1} = g;
for i = 2:4
    k{i} = g*k{i};
    color = colorSpace(90+round(10*activation(i)),:);
    if i==4
        h{i} = drawActiveFrameEllipse(k{i},l_a,er,color);
    else
        h{i} = drawActiveFrameEllipse(k{i},l_b/4,er,color);
    end
end

er=1;
color = colorSpace(90+round(54*activation(n+1)),:);
h{n+1} = drawActiveFrameEllipse(g*gFL,l_a/5,er,color);
color = colorSpace(90+round(54*activation(n+2)),:);
h{n+2} = drawActiveFrameEllipse(g*gFR,l_a/5,er,color);
color = colorSpace(90+round(54*activation(n+3)),:);
h{n+3} = drawActiveFrameEllipse(g*gHL,l_a/5,er,color);
color = colorSpace(90+round(54*activation(n+4)),:);
h{n+4} = drawActiveFrameEllipse(g*gHR,l_a/5,er,color);

g_leg{1}=g*gFL; % should be FR
g_leg{2}=g*gFR; % should be FL
g_leg{3}=g*gHL; % shoule be HR
g_leg{4}=g*gHR; % should be HL
save('g_leg_initial','g_leg');

end

