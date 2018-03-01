iteration=1;
lfs=0.25;
interPhi=pi;
meterPhi=lfs*2*pi;
outerPhi=pi+meterPhi;
if iteration==1
    dutyf=0.75;
    alpha1_col=linspace(0,2*pi,21);
    alpha2_col=linspace(0,2*pi,21);
    g_x_col=zeros(20,20);
    g_y_col=zeros(20,20);
    g_theta_col=zeros(20,20);
    for ind1=19
        for ind2=19
            alpha=@(t) [-pi/6*sin(t+alpha1_col(ind1));-pi/6*sin(t+alpha1_col(ind2))];
            new_drawsalamander
                g_x_col(ind1,ind2)=g_x;
                g_y_col(ind1,ind2)=g_y;
                g_theta_col(ind1,ind2)=g_theta;
        end
    end
elseif iteration==2
    alpha1_col=linspace(0,2*pi,20);
    alpha2_col=linspace(0,2*pi,20);    
    alpha=@(t) [-pi/12*sin(t+alpha1_col(17));-pi/12*sin(t+alpha1_col(17))];
    dutyf_col=0.5:0.05:0.95;
    for dutyf_ind=8
        dutyf=dutyf_col(dutyf_ind);
        new_drawsalamander
        g_x_col(dutyf_ind)=g_x;
        g_y_col(dutyf_ind)=g_y;
        g_theta_col(dutyf_ind)=g_theta;
    end
    plot(dutyf_col,(g_x.^2+g_y.^2).^0.5);
end
    