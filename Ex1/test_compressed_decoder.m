clc % clear command window
clear % clear workspace
close all % close all open plots

% Load image
load image.mat;
tx = signal; % full signal


% Add noise
SNR_dB = 30
wx = gen_noise(SNR_dB, tx);
rx = tx + wx; %received signal

% Demap image signal
rxbits(1:2:2*length(rx)) = uint8(sign(real(rx)) + 1)/2;
rxbits(2:2:2*length(rx)) = uint8(sign(imag(rx)) + 1)/2;


% Show image
compressed_decoder(rxbits, image_size)