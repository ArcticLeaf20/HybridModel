function reactivity = master_reactivity_function(Mod_Reactivity_R, Fuel_Reactivity_R,reac_eff_Xe,u,wierd,o,p)

reactivity = control_rod_reac(u,wierd,o,p) + reac_eff_Xe + Mod_Reactivity_R+ Fuel_Reactivity_R;
        
end 