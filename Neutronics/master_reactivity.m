
if time==0 
        NI=0;
        NX=0;
 end

reactivity_output=master_reactivity_function(1,0.5,NI,NX); %[Power in MW, time in seconds]


function [reac_eff_Xe, NI, NX] = xenon_buildup(power, time, NI,NX)
    ids=24120;           % Iodine_decay_constant [s]
    xds=33120;           % Xenon_decay_constant [s]
    ify=0.0629;          % Iodine_fission_yield [%]
    nu=2.3;               % neutrons per fission 
    xmacs=2.75*10^(-18);  % Xenon_microscopic_abs_cross_section m[cm^2]

    mfcs=0.0487;         % macroscopic_fission_cross_section_U235 [cm^2]
    mifcs=585*10^(-24);    % microscopic_fission_cross_section_U235 [cm^2]
    volume_of_core=5.13*10^27; % [atoms/core]
    E_r=200.7*10^6*1.6*10^(-19);                            % Energy per fission (estimate) converted from MeV to W [W/fisison]
    
    
    iodine_pop=NI;   
    xenon_pop=NX;
   

 

    neutron_flux=power*(10^6)/((mifcs)*E_r*volume_of_core);  % 'Mystery number' is a conversion from MW to W result is in neutrons cm^-2 s^-1
    
    NI =((ify)*(mfcs)*(neutron_flux)-(1/ids)*(iodine_pop))*(time);                              % Determines Iodine population

    iodine_pop = NI;

    NX =((1/ids)*(iodine_pop)-(1/xds)*(xenon_pop)-(xmacs)*(xenon_pop)*(neutron_flux))*(time);  % Determines Xenon population

    xenon_pop = NX;

    reac_eff_Xe = (xmacs*xenon_pop)/(mfcs*nu);  % Calcualtes effect on reactivity from Xenon population
    disp(reac_eff_Xe)
end 
   

function cr_reac_eff = control_rod_reac() % CAN BE FURTHER REFINED, roughly dtermines reactivity; tables with values will later be added for better accuracy
    reg_rod=linspace(0,3.19, 1001);
    safety_rod=linspace(0,2.15,1001);
    shim_rod=linspace(0,2.56,1001);
    trans_rod=linspace(0,2.68,1001);
    reg_height=21;
    safety_height=25;
    shim_height=56;
    trans_height=78;
    cr_reac_eff=(reg_rod(reg_height*10)+safety_rod(safety_height*10)+shim_rod(shim_height*10)+trans_rod(trans_height*10))*0.006;
    disp(cr_reac_eff)

end 



 

function reactivity = master_reactivity_function(power, time, NI, NX)
        reactivity=control_rod_reac()-xenon_buildup(power, time, NI, NX);
        disp(reactivity)
end 
