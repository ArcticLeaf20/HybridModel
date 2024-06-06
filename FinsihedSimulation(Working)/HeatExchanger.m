function PrimaryWaterTout = HeatExchanger(TB,Toutold)
   PrimaryWaterTout = -(((((164.3) * (TB) - 4436.2) / (TB - 12.77)) * ((TB + Toutold) / 2)) / 164.3) + TB;  
end