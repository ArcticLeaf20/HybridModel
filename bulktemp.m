function TB = bulktemp(TB, initial_power,Toutold)
    TB = (178.3 * coreout(TB,initial_power,Toutold) + (341.7 * HeatExchanger(TB,Toutold))) / 520;
   
end