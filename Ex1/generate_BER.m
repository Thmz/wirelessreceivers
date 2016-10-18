clc % clear command window
clear % clear workspace
close all % close all open plots

%% For original constellation points of figure 1.6

nbBits = 500000; % number of bits 

SNR_values = -6:0.5:12; % evaluated SNR values
BER_values = []; % will be calucalated

for SNR_dB = SNR_values
    % Source: Generate random bits
    txbits = randi([0 1],nbBits,1);

    % Map it
    tx = ((txbits(1:2:length(txbits)) - 0.5) + (txbits(2:2:length(txbits)) - 0.5)*i ) * sqrt(2);

    % Add noise 
    wx = gen_noise( SNR_dB, tx);
    rx = tx + wx; %received signal   


    % Demapping: Symbols to bits
    rxbits(1:2:2*length(rx), 1) = uint8(sign(real(rx)) + 1)/2;
    rxbits(2:2:2*length(rx), 1) = uint8(sign(imag(rx)) + 1)/2;
   
    
    % BER: Count errors
    cpr = rxbits ~= txbits;
    errors = sum(cpr);
    
    % Add BER to matrix
    ber = errors/nbBits;
    BER_values = [ BER_values ber];
end


%% For alternative constellation points of figure 1.6

BER_values_alternative = []; % will be calucalated
for SNR_dB = SNR_values
    % Source: Generate random bits
    txbits = randi([0 1],nbBits,1);

    % Map it
    tx = exp(1i * (txbits - 0.5 ) * pi/2);
    
    % Add noise 
    wx = gen_noise( SNR_dB, tx);
    rx = tx + wx; %received signal   


    % Demapping: Symbols to bits
    rxbits = mod(floor((angle(rx) + pi/2)/(pi/2)), 4);
    
    % BER: Count errors
    cpr = rxbits ~= txbits;
    errors = sum(cpr);
    
    % Add BER to matrix
    ber = errors/nbBits;
    BER_values_alternative = [ BER_values_alternative ber];
end

%% Make plot
figure('Name', 'BER');
semilogy(SNR_values, BER_values);
hold on;
semilogy(SNR_values, BER_values_alternative);
xlabel('SNR [dB]');
ylabel('BER [%]');
title('BER');
legend('With original mapping','Alternative mapping of figure 1.6');