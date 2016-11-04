% task_1_and_2.m: Example solution for Tasks 1.4.1 and 1.4.2
%
% Author: Alexios Balatsoukas-Stimming
%

clear all
close all
clc

% SNR
SNR = 10;

% Load file
load image

% Convert SNR from dB to linear
SNRlin = 10^(SNR/10);

% Add AWGN
noisy_signal = signal + sqrt(1/(2*SNRlin)) * (randn(size(signal)) + 1i*randn(size(signal)) );  

% Demap
img = demapper(noisy_signal);

% Display image
image_decoder(img,image_size);

% Plot noisy symbols on complex plane
plot(noisy_signal, '.');
title('Noisy received symbols')
xlabel('Real part')
ylabel('Imaginary part')