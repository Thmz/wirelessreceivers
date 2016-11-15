%function ex3_3_1_demodpulse
% ex3_3_1 - Task 3.3.1
%    Demodulate root raised cosine shaped pulses
%    
% 
% Author(s): Nicholas Preyss
% Copyright (c) 2012 TCL.
% 

SNR           = 10;
rx_filterlen  = 12; % length of receive filter
os_factor     = 4;  % oversampling factor

% load shaped symbols
load ./task3.mat
data_length = prod(image_size) * 8 / 2;

Pa = sum(abs(tx).^2)/length(tx) % calculate signal power

% convert SNR from dB to linear
SNRlin = 10^(SNR/10);

% add awgn channel
wx =  1/sqrt(SNRlin*2) * (randn(length(signal), 1) + randn(length(signal), 1)*1i );
rx_signal = signal + wx;

% apply matched filter
fir_rrc =  rrc(os_factor, 0.22, rx_filterlen);
filtered_rx_signal = conv(rx_signal.', fir_rrc, 'same');

% find start of data frame
beginning_of_data = frame_sync(filtered_rx_signal.', os_factor); % Index of the first data symbol

% decode image
image_decoder(demapper(filtered_rx_signal(beginning_of_data : os_factor : beginning_of_data + os_factor * (data_length - 1))), image_size);