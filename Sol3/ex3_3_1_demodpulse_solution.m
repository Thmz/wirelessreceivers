function ex3_3_1_demodpulse_solution
% ex3_3_1 - Task 3.3.1 solution
%    Demodulate root raised cosine shaped pulses
%    
% 
% Author(s): Nicholas Preyss
% Copyright (c) 2012 TCL.
% 

SNR           = 10;
rx_filter     = 12; % length of receive filter
os_factor     = 4;  % oversampling factor

% load shaped symbols
load ./task3.mat
data_length = prod(image_size) * 8 / 2;

% convert SNR from dB to linear
SNRlin = 10^(SNR/10);

% add awgn channel
rx_signal = signal + sqrt(1/(2*SNRlin)) * (randn(size(signal)) + 1i*randn(size(signal)) ); 

% apply matched filter
rolloff = 0.22;
pulse               = rrc(os_factor,rolloff,rx_filter);
filtered_rx_signal  = zeros(1,length(rx_signal)-2*rx_filter);
data                = rx_signal(:);

for i=1:length(data)-2*rx_filter
    
    segment   = data(i:i+2*rx_filter);
    filtered_rx_signal(i) = segment.'*pulse;
    
end

% find start of data frame
beginning_of_data = frame_sync(filtered_rx_signal.', os_factor); % Index of the first data symbol

% decode image
image_decoder(demapper(filtered_rx_signal(beginning_of_data : os_factor : beginning_of_data + os_factor * (data_length - 1))), image_size);