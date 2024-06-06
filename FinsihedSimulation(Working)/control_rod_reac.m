function cr_reac_eff = control_rod_reac(u,wierd,o,p) % CAN BE FURTHER REFINED, roughly dtermines reactivity; tables with values will later be added for better accuracy
    reg_rod=linspace(0,3.19, 1001);
    safety_rod=linspace(0,2.15,1001);
    shim_rod=linspace(0,2.56,1001);
    trans_rod=linspace(0,2.68,1001);
    if u>100
        u=100;
    end
    if wierd>100
        wierd=100;
    end
    if o>100
        o=100;
    end
    if p>100
        p=100;
    end
    reg_height=p;
    safety_height=wierd;
    shim_height=o;
    trans_height=u;
   
    %cr_reac_eff = -1*(old_Mod_Reactivity_R+old_Fuel_Reactivity_R);
    cr_reac_eff=(reg_rod(1001-reg_height*10)+safety_rod(1001-(safety_height*10))+shim_rod(1001-shim_height*10)+trans_rod(1001-trans_height*10))*0.0075;

end 