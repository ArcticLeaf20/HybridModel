function modtemp = avgmodtemp(TB, initial_power,Toutold)
    modtemp = (HeatExchanger(TB,Toutold) + coreout(TB, initial_power,Toutold)) / 2;
end