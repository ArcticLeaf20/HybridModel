function cr_reac_eff = control_rod_reac(old_Mod_Reactivity_R,old_Fuel_Reactivity_R) % CAN BE FURTHER REFINED, roughly dtermines reactivity; tables with values will later be added for better accuracy
    reg_rod=linspace(0,3.19, 1001);
    safety_rod=linspace(0,2.15,1001);
    shim_rod=linspace(0,2.56,1001);
    trans_rod=linspace(0,2.68,1001);
    reg_height=89;
    safety_height=91;
    shim_height=94;
    trans_height=98;
    %cr_reac_eff = -1*(old_Mod_Reactivity_R+old_Fuel_Reactivity_R);
    cr_reac_eff=(reg_rod(1000-reg_height*10)+safety_rod(1000-safety_height*10)+shim_rod(1000-shim_height*10)+trans_rod(1000-trans_height*10))*0.0075;

end 