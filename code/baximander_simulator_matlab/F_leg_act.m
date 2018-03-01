function leg_act_col=F_leg_act(time,dutyf)
leg_act_col=zeros(size(time));
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
        leg_act=0;
    elseif t<2*pi-(1-dutyf)*pi
        leg_act=1;
    else
        leg_act=0;
    end
    leg_act_col(t_ind)=leg_act;
end

