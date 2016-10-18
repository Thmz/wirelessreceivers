clc % clear command window
clear % clear workspace
close all % close all open plots

% Loop for multiple SNR's
SNR_values = [3 6 20 60]

for SNR_dB = SNR_values
    
    % Load image
    load image.mat;
    tx = signal; % full signal 
    
    
    % Add noise
    wx = gen_noise(SNR_dB, tx);
    rx = tx + wx; %received signal
    snr(rx, wx)
    
    % Plot constellation points
    figure('Name', strcat('Constellation points for SNR ', num2str(SNR_dB), ' dB'));
    plot(rx(1:10000), '.'); % limit length of signal for constellation points plot
    title(strcat('Constellation points for SNR ', num2str(SNR_dB), ' dB'));
    xlabel('Real');
    ylabel('Imaginary');
    
    % Demap image signal
    rxbits(1:2:2*length(rx)) = uint8(sign(real(rx)) + 1)/2;
    rxbits(2:2:2*length(rx)) = uint8(sign(imag(rx)) + 1)/2;
    
    
    % Show original image
    image_decoder(rxbits, image_size);
    %compressed_decoder(rxbits, image_size)
    
    
end
