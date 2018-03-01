function alpha3=get_apha(t)
alpha3=0;
if t<pi/4 
    alpha3=pi/3*sin(t*2);
end
if t>pi/4 && t<pi*2/4
    alpha3=pi/3*cos((t-pi/4)*2);
end
if t>pi*2/4 && t<pi*3/4
    alpha3=3*pi/4*cos((t-pi/4)*2);
end
if t>pi*3/4 && t<pi*4/4
    alpha3=-3*pi/4;
end
if t>pi*4/4 && t<pi*5/4
    alpha3=-3*pi/4+pi/4*sin((t-pi)*2);
end
if t>pi*5/4 && t<pi*7/4
    alpha3=-pi/2+pi/2*sin((t-pi*5/4));
end