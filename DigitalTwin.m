format default



time_values_real = linspace(1, 604800,604800);                   % Initializes Primary Time loop (wraps time around each loop)

PRKEtime_values = linspace(0, 0.001, 10000);
generationTime = 1 * (10^-7);                                                   % Mean neutron generation time
decay_constants = [0.0128, 0.0301, 0.124, 0.325, 1.12, 2.69];                   % Decay constants for each precursor group
beta_i = [0.000073, 0.000626, 0.000443, 0.000685, 0.000181, 0.000092];          % Beta values for precursor groups
beta = sum(beta_i);                                                             % Effective beta 
concentration_values = zeros(6, length(PRKEtime_values));                       % Array for concentration values
TB=34;                                                                          % Initial bulk temperature
initial_power = 1;                                                              % Initial power in MW
Toutold=27;                                                                     % Initial Tmeperature out of HX
reactivity=0;                                                                   % Initial Reactivity 
old_Mod_Reactivity_R = -0.0015;                                                 % Initial Moderator effect on reactivity

old_Fuel_Reactivity_R = -0.0038;                                                % Initial Fuel Temp effect on reactivity

initial_concentration = beta_i./(decay_constants.*generationTime).*initial_power;         % Initializes concentration values


ids=24120;            % Iodine_decay_constant [s]                               
xds=33120;            % Xenon_decay_constant [s]
ify=0.0629;           % Iodine_fission_yield [%]


nu=2.3;               % neutrons per fission 
xmacs=2.75*10^(-18);  % Xenon_microscopic_abs_cross_section m[cm^2]

mfcs=0.0487;         % macroscopic_fission_cross_section_U235 [cm^2]
mifcs=585*10^(-24);    % microscopic_fission_cross_section_U235 [cm^2]
volume_of_core=5.13*10^27;                                                      % [atoms/core]
E_r=200.7*10^6*1.6*10^(-19);                                                    % Energy per fission (estimate) converted from MeV to W [W/fisison]

time_values = linspace(0, 0.01, 40);
reac_values = zeros(1, length(time_values));

initial_iodine_pop=0;
initial_xenon_pop=0;


dt =PRKEtime_values(2)-PRKEtime_values(1);
concentration_sum = 0;

actual_time=[];
tot_p=[];
TB_tot=[];
f_r_grap=[];
for x = time_values_real


for i = 1:length(PRKEtime_values)                % Euler method that updates power value each 0.001 second
    
    concentration_values(:,i) = initial_concentration;

    concentration_sum = decay_constants*initial_concentration';
    A = (reactivity - beta) / generationTime;
   
    k1 = initial_power*A+concentration_sum;
    % k2 = (initial_power + dt*(k1/2)) * (exp(generationTime*(dt/2))) + (1/A) * (exp(generationTime*(dt/2))-1) * concentration_sum;
    % k3 = (initial_power + dt*(k2/2)) * (exp(generationTime*(dt/2))) + (1/A) * (exp(generationTime*(dt/2))-1) * concentration_sum;
    % k4 = (initial_power + dt*k3) * (exp(generationTime*dt)) + (1/A) * (exp(generationTime*dt)-1) * concentration_sum;
    % new_power = (initial_power + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)); %RK4
    new_power = initial_power+dt*k1;

    k1_c = beta_i./generationTime.*initial_power-decay_constants.*initial_concentration;
    new_concentration = initial_concentration+dt*k1_c;
    % for z = 1:6
    % concentration_values(z) = initial_concentration(z) * exp(-1*decay_constants(z)*dt) + (beta / generationTime) * (1/2)*(initial_power+old_p) * (1-exp(-1*decay_constants(z)*dt)) * (1/decay_constants(z));
    % end
    initial_power = new_power;
    initial_concentration = new_concentration;
    
    

end






%updates moderator temp in core with new T bulk value
Tf = Fueltemperature(TB, initial_power,Toutold);

modtemp = avgmodtemp(TB, initial_power,Toutold);

Mod_Reactivity_R = mod_reactivity(modtemp);

Fuel_Reactivity_R = Fuel_reactivity(Tf);

Toutold =HeatExchanger(TB, Toutold);   %updates temp with new power

TB = bulktemp(TB, initial_power,Toutold);   %updates temp with new power
for i = 1:length(PRKEtime_values)  % Determines Xenon population as time progresse also using Euler
    
    
    neutron_flux = (initial_power*(10^6))/((mifcs)*E_r*volume_of_core);

    
    k1_i =((ify) * (mfcs)*(neutron_flux)-(1/ids)*(initial_iodine_pop));                              % Determines Iodine population
    
    new_iodine= initial_iodine_pop + dt*k1_i;
    
    
    k1_x =((1/ids)*(initial_iodine_pop)-(1/xds)*(initial_xenon_pop)-(xmacs)*(initial_xenon_pop)*(neutron_flux));  % Determines Xenon population

    new_xenon = initial_xenon_pop + dt*k1_x;

    
    initial_iodine_pop = new_iodine;
    initial_xenon_pop = new_xenon;

 
    
end 

reac_eff_Xe = (xmacs*new_xenon)/(mfcs*nu);  % Calculates effect on reactivity from Xenon population


reactivity = master_reactivity_function(old_Mod_Reactivity_R, old_Fuel_Reactivity_R, Mod_Reactivity_R, Fuel_Reactivity_R, reac_eff_Xe);     % Determines overall rectivity (updates reactivty for next PRKE iteration)

old_Mod_Reactivity_R = Mod_Reactivity_R;

old_Fuel_Reactivity_R = Fuel_Reactivity_R;    % Fudge values for control rod height for now to make system not converge


% Everything left in loop is for graphing the power or other paremters as need be
actual_time(x) = x/1000;    

tot_p(x)=initial_power;

%TB_tot(x)=TB;

%f_r_graph(x)=Fuel_Reactivity_R


%plot(actual_time,f_r_grap),grid,
%xlabel('Time (s)'), ylabel('Fuel Reactivty')
%plot(actual_time,TB_tot),grid,
%xlabel('Time (s)'), ylabel('Temp (C)')


plot(actual_time,tot_p),grid,
xlabel('Time (s)'), ylabel('Power (MW)')
drawnow

disp(initial_power)


end 
   




















% Fuctions used in loops above, should be self explanatory


function Mod_Reactivity_R = mod_reactivity(T_avg_in_core)
T = T_avg_in_core; % in celcius
B_eff = 0.0075;
MTC = -0.0072; % dollars per celcius
MR_D = MTC*T; % moderator reactivity in cents
Mod_Reactivity_R = MR_D*B_eff;
end 



function Fuel_Reactivity_R = Fuel_reactivity(T_avg_fuel_in_core)
Fuel_Reactivity_R = -1.072*10^-5*(T_avg_fuel_in_core)+2.49*10^-3; % in celcius
end 







function PrimaryWaterTout = HeatExchanger(TB,Toutold)
   PrimaryWaterTout = -(((((164.3) * (TB) - 4436.2) / (TB - 12.77)) * ((TB + Toutold) / 2)) / 164.3) + TB;  
end

function Tout = coreout(TB, initial_power, Toutold)
    Tout = HeatExchanger(TB,Toutold) + ((initial_power*1000) / (39.12 * 4.12));
end


function Tf = Fueltemperature(TB, initial_power,Toutold)
    Tf = 558.895 + coreout(TB, initial_power,Toutold);
    
end

function modtemp = avgmodtemp(TB, initial_power,Toutold)
    modtemp = (HeatExchanger(TB,Toutold) + coreout(TB, initial_power,Toutold)) / 2;
end


function TB = bulktemp(TB, initial_power,Toutold)
    TB = (178.3 * coreout(TB,initial_power,Toutold) + (341.7 * HeatExchanger(TB,Toutold))) / 520;
   
end

function cr_reac_eff = control_rod_reac(old_Mod_Reactivity_R,old_Fuel_Reactivity_R) % CAN BE FURTHER REFINED, roughly dtermines reactivity; tables with values will later be added for better accuracy
    reg_rod=linspace(0,3.19, 1001);
    safety_rod=linspace(0,2.15,1001);
    shim_rod=linspace(0,2.56,1001);
    trans_rod=linspace(0,2.68,1001);
    reg_height=36;
    safety_height=5;
    shim_height=5;
    trans_height=6;
    %cr_reac_eff=(reg_rod(reg_height*10)+safety_rod(safety_height*10)+shim_rod(shim_height*10)+trans_rod(trans_height*10))*0.0075;
    cr_reac_eff = -1*(old_Mod_Reactivity_R+old_Fuel_Reactivity_R);
    

end 




function reactivity = master_reactivity_function(old_Mod_Reactivity_R,old_Fuel_Reactivity_R,Mod_Reactivity_R, Fuel_Reactivity_R,reac_eff_Xe)
        reactivity = control_rod_reac(old_Mod_Reactivity_R,old_Fuel_Reactivity_R) + reac_eff_Xe*0 + Mod_Reactivity_R+ Fuel_Reactivity_R
        
end 


