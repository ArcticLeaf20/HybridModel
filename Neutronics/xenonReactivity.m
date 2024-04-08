% OSTR Hybrid model, Xenon reactivity module
% Calculates effect of xenon on reactivity for use in PRKE

% Parameters

ids = 24120;           % Iodine_decay_constant [s]
xds = 33120;           % Xenon_decay_constant [s]
ify = 0.0629;          % Iodine_fission_yield [%]
   
xmacs = 3.5*10^(-18);  % Xenon_microscopic_abs_cross_section m[cm^2]
flux = neutron_flux;
mfcs = 0.0487;         % macroscopic_fission_cross_section [cm^2]
mifcs= 1*10^(-24);     % microscopic_fission_cross_section [cm-1]
t = time;


NI=0;                             
iodine_pop=NI;
NX=0;   
xenon_pop=NX;


xenon_buildup [change_in_reactivity] = flux(neutron_flux, time)
     NI =((ify)*(mfcs)*(flux)-(1/ids)*(iodine_pop))*(time);                              % Determines Iodine population
     iodine_pop = NI;
     NX =((1/ids)*(iodine_pop)-(1/xds)*(xenon_pop)-(xmacs)*(xenon_pop)*(neutron_flux))*(time);   % Determines Xenon population
     xenon_pop = NX;


     reac_eff = ((xmacs)*(xenon_pop*neutron_flux*time))/((mifcs)*(neutron_flux*time));  % Calcualtes effect on reactovty from Xenon population
