%clear all
%close all
%clc

% Set parameters
os_factor = 4;
SNR = 30;

load ber_pn_seq
load pn_sequence

data_length = length(ber_pn_seq)/2; % Number of QPSK data symbols

% Matched filter
filtered_rx_signal = matched_filter(awgn(signal, SNR), os_factor, 6); % 6 is a good value for the one-sided RRC length (i.e. the filter has 13 taps in total)

% Frame synchronization
data_idx = frame_sync(filtered_rx_signal, os_factor); % Index of the first data symbol

data = zeros(1,data_length);
data2 = zeros(1,data_length);

cum_err = 0;
diff_err = zeros(1,data_length);
epsilon  = zeros(1,data_length);

frame_sync_length = 100;

% Use preamble symbols to improve timing offset estimation (Task 3)
for i = floor(data_idx/os_factor)-frame_sync_length:floor(data_idx/os_factor)-1
    
    idx_start  = i*os_factor;
     
    idx_range  = idx_start:idx_start+os_factor-1;
    segment    = filtered_rx_signal(idx_range);
    
    % Estimate timing error epsilon
    pwr         = abs(segment).^2;
    diff_err = [1 -1j -1 1j]*pwr;
    cum_err     = cum_err + diff_err;
    
end

for i=1:data_length
    
     idx_start  = data_idx+(i-1)*os_factor;
     
     idx_range  = idx_start:idx_start+os_factor-1;
     segment    = filtered_rx_signal(idx_range);
    
     % Estimate timing error epsilon
     pwr         = abs(segment).^2;
     diff_err(i) = [1 -1j -1 1j]*pwr;
     cum_err     = cum_err + diff_err(i);
     epsilon(i)  = -1/(2*pi)*angle(cum_err);
     
     % Interpolate
     sample_diff   = floor(epsilon(i)*os_factor); % integer
     int_diff      = mod(epsilon(i)*os_factor,1); % interval [0 1)
    
     
     % linear
     y     = filtered_rx_signal(idx_start+sample_diff:idx_start+sample_diff+1);
     y_hat = y(1)+int_diff*(y(2)-y(1));
     data(i) = y_hat;
     
     % cubic
     y2     = filtered_rx_signal(idx_start+sample_diff-1:idx_start+sample_diff+2);
     y_hat2 = cubic_interpolator(y2,int_diff);
     data2(i) = y_hat2;     

     
end

% Plot epsilon
plot(1:data_length, epsilon)

BER_lin = mean(demapper(data) ~= ber_pn_seq)
BER_cub = mean(demapper(data2) ~= ber_pn_seq)
