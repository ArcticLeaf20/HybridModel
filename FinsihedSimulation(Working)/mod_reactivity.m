function Mod_Reactivity_R = mod_reactivity(T_avg_in_core)
T = T_avg_in_core; % in celcius
B_eff = 0.0075;
MTC = -0.0072; % dollars per celcius -0.0072 to -0.004
MR_D = MTC*T; % moderator reactivity in cents
Mod_Reactivity_R = MR_D*B_eff;
end 