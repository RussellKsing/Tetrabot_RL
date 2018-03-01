function xi=questions(M1,N1,M2,N2,M3,N3,xi_pre)
theta_o=xi_pre(1);
a1=M1(1,1);b1=M1(2,1);c1=M1(1,3);d1=M1(2,3);e1=N1(1,1);f1=N1(2,1);
a2=M2(1,1);b2=M2(2,1);c2=M2(1,3);d2=M2(2,3);e2=N2(1,1);f2=N2(2,1);
a3=M3(1,1);b3=M3(2,1);c3=M3(1,3);d3=M3(2,3);e3=N3(1,1);f3=N3(2,1);


g1=(a1*b2-b1*a2);h1=(c1*b2-d1*a2-c2*b2+d2*a2);i1=b2*(e1-e2)-a2*(f1-f2);j1=b2*(f1-f2)-a2*(e2-e1);
g2=(a1*b3-b1*a3);h2=(c1*b3-d1*a3-c3*b3+d3*a3);i2=b3*(e1-e3)-a3*(f1-f3);j2=b3*(f1-f3)-a3*(e3-e1);

myfun=@(theta) (i1*g2-i2*g1)*cos(theta)+(j1*g2-j2*g1)*sin(theta)-(h1*g2-h2*g1);
theta= fsolve(@(x) myfun(x),theta_o);

z_r=[(e1-e2)*cos(theta)+(f1-f2)*sin(theta)-c1+c2;
     (e1-e3)*cos(theta)+(f1-f3)*sin(theta)-c1+c3;
     (e2-e1)*sin(theta)+(f1-f2)*cos(theta)-d1+d2;];
z_r(2)=z_r(1)+z_r(2);
z_l=[a1    -a2   0;
     2*a1  -a2 -a3;
     b1    -b2   0;];
Z=inv(z_l)*z_r;
z1=Z(1);z2=Z(2);z3=Z(3);

x_l=[e1*cos(theta)+f1*sin(theta)-z1*a1-c1;
    -e1*sin(theta)+f1*cos(theta)-z1*b1-d1;];
x_r=[cos(theta) sin(theta)
    -sin(theta) cos(theta)];
X=inv(x_r)*x_l;
x=X(1);y=X(2);
xi=[theta x y z1 z2 z3];
end
        