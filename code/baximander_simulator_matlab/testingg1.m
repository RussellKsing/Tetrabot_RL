
clc;clear;
dutyf = 0.9;
interPhi = 4.01212602265;
meterPhi = 1.18943521037;
outerPhi = -1.41714577328;
a2_1 = 0.701045249749;
a2_2 = 0.905461661386;
b2_1 = 2.86276125191;
b2_2 = 0.0474152993261;
c_FL = 0.900959302491;
c_FR = 0.790483551065;
c_HL = 0.84577315574;
c_HR = 0.606085753847;


l_a = 2.25;
l_b = 12.1;
l_c = 4.55;
l_d = 6.8;
l_e = 10.35;


for phase1_ind=17
    for phase2_ind=17
        %phase1=phase1_col(phase1_ind);
        %phase2=phase2_col(phase2_ind);
        alpha=@(t) [a2_1*sin(t+b2_1);a2_2*sin(t+b2_2);0];
        %alpha=@(t) [0;0;get_apha(t)];
        new_drawsalamander
    end
end

