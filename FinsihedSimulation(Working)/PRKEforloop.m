
function [initial_power,reactivity,init1,init2,init3,init4,init5,init6,Toutold,TB,initial_iodine_pop,initial_xenon_pop,Tf] = PRKEforloop(initial_power,reactivity,init1,init2,init3,init4,init5,init6,Toutold,TB,initial_iodine_pop,initial_xenon_pop,u,wierd,o,p)
format long
PRKEtime_values = linspace(0, 0.01, 10000);
generationTime = 1 * (10^-7); % Mean neutron generation time
decay_constants = [0.00124, 0.031, 0.111, 0.301, 1.136, 3.014]; % Decay constants for each precursor group
beta_i = [0.00021, 0.00141, 0.00127, 0.00255, 0.00074, 0.00027]; % Beta values for precursor groups
beta = sum(beta_i); % Effective beta;
concentration_values = zeros(6, length(PRKEtime_values)); % Array for concentration values
initial_concentration=[init1,init2,init3,init4,init5,init6];





ids=24120;            % Iodine_decay_constant [s]
xds=33120;            % Xenon_decay_constant [s]
ify=0.0629;           % Iodine_fission_yield [%]


nu=2.3;               % neutrons per fission 
xmacs=2.75*10^(-18);  % Xenon_microscopic_abs_cross_section m[cm^2]

mfcs=0.0487;         % macroscopic_fission_cross_section_U235 [cm^2]
mifcs=585*10^(-24);    % microscopic_fission_cross_section_U235 [cm^2]
volume_of_core=5.13*10^27; % [atoms/core]
E_r=200.7*10^6*1.6*10^(-19);                            % Energy per fission (estimate) converted from MeV to W [W/fisison]

time_values = linspace(0, 0.01, 40);
reac_values = zeros(1, length(time_values));




dt =PRKEtime_values(2)-PRKEtime_values(1);


concentration_sum=0;

for i = 1:length(PRKEtime_values)
    
    concentration_values(:,i) = initial_concentration;

    concentration_sum = decay_constants*initial_concentration';
    A = (reactivity - beta) / generationTime;
   
    k1 = initial_power*A+concentration_sum;
    k2 = (initial_power + dt*(k1/2)) *A+concentration_sum;
    k3 = (initial_power + dt*(k2/2)) *A+concentration_sum;
    k4 = (initial_power + dt*k3) *A+concentration_sum;
    % new_power = (initial_power + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)); %RK4
    new_power = initial_power + (dt/6) * (k1+k2+k3+k4);

    k1_c = beta_i./generationTime.*initial_power-decay_constants.*(initial_concentration);
    k2_c = beta_i./generationTime.*initial_power-decay_constants.*(initial_concentration + dt*(k1_c/2));
    k3_c = beta_i./generationTime.*initial_power-decay_constants.*(initial_concentration+ dt*(k2_c/2));
    k4_c = beta_i./generationTime.*initial_power-decay_constants.*(initial_concentration+ dt*k3_c);
    new_concentration = initial_concentration+dt+ (dt/6) * (k1_c+k2_c+k3_c+k4_c);
    % for z = 1:6
    % concentration_values(z) = initial_concentration(z) * exp(-1*decay_constants(z)*dt) + (beta / generationTime) * (1/2)*(initial_power+old_p) * (1-exp(-1*decay_constants(z)*dt)) * (1/decay_constants(z));
    % end
    initial_power = new_power;
    initial_concentration = new_concentration;
    
    

end

init1=initial_concentration(1);
init2=initial_concentration(2);
init3=initial_concentration(3);
init4=initial_concentration(4);
init5=initial_concentration(5);
init6=initial_concentration(6);

%updates moderator temp in core with new T bulk value
Tf = Fueltemperature(TB, initial_power,Toutold);

modtemp = avgmodtemp(TB, initial_power,Toutold);

Mod_Reactivity_R = mod_reactivity(modtemp);

Fuel_Reactivity_R = Fuel_reactivity(Tf);

Toutold =HeatExchanger(TB, Toutold);

TB = bulktemp(TB, initial_power,Toutold)+(-1+(1+1)*rand(1,1));
for i = 1:length(PRKEtime_values)
    
    
    neutron_flux = (initial_power*(10^6))/((mifcs)*E_r*volume_of_core);

    
    k1_i =((ify) * (mfcs)*(neutron_flux)-(1/ids)*(initial_iodine_pop));                              % Determines Iodine population
    
    new_iodine= initial_iodine_pop + dt*k1_i;
    
    
    k1_x =((1/ids)*(initial_iodine_pop)-(1/xds)*(initial_xenon_pop)-(xmacs)*(initial_xenon_pop)*(neutron_flux));  % Determines Xenon population

    new_xenon = initial_xenon_pop + dt*k1_x;

    
    initial_iodine_pop = new_iodine;
    initial_xenon_pop = new_xenon;

 
    
end 

reac_eff_Xe = (xmacs*new_xenon)/(mfcs*nu);  % Calculates effect on reactivity from Xenon population


reactivity = master_reactivity_function(Mod_Reactivity_R, Fuel_Reactivity_R, reac_eff_Xe,u,wierd,o,p);




end 