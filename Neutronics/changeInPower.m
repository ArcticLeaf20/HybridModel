function dPdt = changeInPower(time, power, rho, concentration)
    dPdt = ( (rho - beta) / generationTime) * power(time) + (lambda * concentration(time));
end

