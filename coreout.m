function Tout = coreout(TB, initial_power, Toutold)
    Tout = HeatExchanger(TB,Toutold) + ((initial_power*1000) / (39.12 * 4.12));
end