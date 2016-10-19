clear 
clc

os_factor = 4;
SNR = 10;

load ./task4

data_length = prod(image_size) * 8 / 2; % Number of QPSK data symbols

% convert SNR from dB to linear
SNRlin = 10^(SNR/10);

% add awgn channel
rx_signal = signal + sqrt(1/(2*SNRlin)) * (randn(size(signal)) + 1i*randn(size(signal)) ); 

% Matched filter
filtered_rx_signal = matched_filter(rx_signal, os_factor, 6); % 6 is a good value for the one-sided RRC length (i.e. the filter has 13 taps in total)

% Frame synchronization
data_idx = frame_sync(filtered_rx_signal, os_factor); % Index of the first data symbol



data_linear = zeros(1,data_length);

cum_err = 0;
diff_err = zeros(1,data_length);
epsilon  = zeros(1,data_length);

for i=1:data_length
    
     idx_start  = data_idx+(i-1)*os_factor;
     idx_range  = idx_start:idx_start+os_factor-1;
     segment    = filtered_rx_signal(idx_range);
     
     
     % Estimate timing error epsilon
     for k =  idx_start : 1 : idx_start + os_factor -1
          cum_err = cum_err + filtered_rx_signal(k) * (-1i)^k;
     end
     epsilon(i) = -1/(2*pi) * angle(cum_err);
     
     
     % Chose best sampling point
     eps_floored = floor(epsilon(i));
     idx_sample = idx_start + eps_floored;
     interval_fraction = epsilon(i) - eps_floored;
     
     % Interpolate     
     data_linear(i) = linear_interpolator(filtered_rx_signal, idx_sample, interval_fraction);
     
end

plot(epsilon)
image_decoder(demapper(data_linear), image_size);