function [ h ] = drawActiveFrameEllipse( g, L,er, color )

beta=linspace(0,2*pi,101);
e_x=L*cos(beta);
e_y=er*L*sin(beta);
e_p=[e_x;e_y];
R=g(1:2,1:2);
e_p=R*e_p;

h = patch(g(1,3)+e_p(1,:),g(2,3)+e_p(2,:),color,'linewidth',1);

end

