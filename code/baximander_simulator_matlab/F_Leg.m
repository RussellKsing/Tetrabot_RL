function leg_angle_col=F_Leg(Aleg,time,dutyf)
leg_angle_col=zeros(size(time));
for t_ind=1:length(time);
    t=time(t_ind);
    while t<0 || t>2*pi
        if t<0
            t=t+2*pi;
        elseif t>2*pi
            t=t-2*pi;
        end
    end
    
    if t<(1-dutyf)*pi;
        leg_angle=Aleg*sin(t/2/(1-dutyf));
    elseif t<2*pi-(1-dutyf)*pi
        leg_angle=Aleg*cos((t-(1-dutyf)*pi)*(pi/(2*pi-2*(1-dutyf)*pi)));
    else
        leg_angle=-Aleg*cos((t-2*pi+(1-dutyf)*pi)/2/(1-dutyf));
    end
    leg_angle_col(t_ind)=leg_angle;
end

