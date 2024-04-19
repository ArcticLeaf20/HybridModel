
reactivity=master_reactivity_function(1,0.5);


function reac_eff_Xe  = xenon_buildup(power, time)
    ids=24120;           % Iodine_decay_constant [s]
    xds=33120;           % Xenon_decay_constant [s]
    ify=0.0629;          % Iodine_fission_yield [%]
   
    xmacs=3.5*10^(-18);  % Xenon_microscopic_abs_cross_section m[cm^2]

    mfcs=0.0487;         % macroscopic_fission_cross_section_U235 [cm^2]
    mifcs=585*10^(-24);    % microscopic_fission_cross_section_U235 [cm^2]
    volume_of_core=5.13*10^27; % [atoms/core]
    E_r=200.7;                            % Energy per fission (estimate) [MeV/fisison]
    h=0.5;

    NI=0;                             
    iodine_pop=NI;
    NX=0;   
    xenon_pop=NX;
    neutron_flux=power/((mifcs)*E_r*volume_of_core*10^6*1.6*10^-19);
    
    NI =((ify)*(mfcs)*(neutron_flux)-(1/ids)*(iodine_pop))*(time);                              % Determines Iodine population

    iodine_pop = NI;

    NX =((1/ids)*(iodine_pop)-(1/xds)*(xenon_pop)-(xmacs)*(xenon_pop)*(neutron_flux))*(time);   % Determines Xenon population

    xenon_pop = NX;

    reac_eff_Xe = ((xmacs)*(xenon_pop*neutron_flux*time))/((mifcs)*(neutron_flux*time));  % Calcualtes effect on reactivity from Xenon population
end 
   

function cr_reac_eff = control_rod_reac(power) % NOT COMPLETE
        cr_reac_eff=power*0;
end 

function reactivity = master_reactivity_function(power, time)
        reactivity=xenon_buildup(power, time)+control_rod_reac(power);
        disp(reactivity)
end 


    
    