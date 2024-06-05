function reactivity = master_reactivity_function(old_Mod_Reactivity_R,old_Fuel_Reactivity_R,Mod_Reactivity_R, Fuel_Reactivity_R,reac_eff_Xe)
        reactivity = control_rod_reac(old_Mod_Reactivity_R,old_Fuel_Reactivity_R) + reac_eff_Xe + Mod_Reactivity_R+ Fuel_Reactivity_R
        
end 