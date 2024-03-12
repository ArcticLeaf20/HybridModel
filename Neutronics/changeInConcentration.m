function dCdt = changeInConcentration(time, power, concentration)
    dCdt = (beta - generationTime) * power(time) - lambda * concentration(time);
end
