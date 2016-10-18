function [ BER ] = calc_ber( txbits, rxbits )
%CALC_BER calculates bit error rate

cpr = rxbits ~= txbits;
errors = sum(cpr);
BER = errors/length(txbits);

end

