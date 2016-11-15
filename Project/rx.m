function [rxbits conf] = rx(rxsignal,conf,k)
% Digital Receiver
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete causal
%   receiver in digital domain.
%
%   rxsignal    : received signal
%   conf        : configuration structure
%   k           : frame index
%
%   Outputs
%
%   rxbits      : received bits
%   conf        : configuration structure
%

% dummy 
%rxbits = zeros(conf.nbits,1);

% Downconvert
time = 1:1/conf.f_s:1+(length(rxsignal)-1)/conf.f_s;
length(rxsignal)
length(time)
rxsignal = rxsignal .* exp(-1j*2*pi*conf.f_c * time.');
rxsignal = lowpass(rxsignal, conf);


data_length = length(rxsignal);


% Interpolator type
interpolator_type = 'linear';

% Receiver
filtered_rx_signal = matched_filter(rxsignal, conf.os_factor, 6); % 6 is a good value for the one-sided RRC length (i.e. the filter has 13 taps in total)


% Frame synchronization
[data_idx theta] = frame_sync(filtered_rx_signal, conf.os_factor) % Index of the first data symbol
data_idxInit = data_idx;

payload_data = zeros(data_length, 1); % The payload data symbols
eps_hat = zeros(data_length, 1);     % Estimate of the timing error
theta_hat = zeros(data_length, 1);   % Estimate of the carrier phase
theta_hat(1) = theta;

% Loop over the data symbols with estimation and correction of the timing
% error and phase

for k = 1 : data_length,
    
    % timing estimator
    %eps_hat(k) = timing_estimator(filtered_rx_signal(data_idx : data_idx + conf.os_factor - 1));
%      for jj =  data_idx : data_idx + conf.os_factor - 1
%           cum_err = cum_err + filtered_rx_signal(jj) * (-1i)^jj;
%      end
%      eps_hat(k) = -1/(2*pi) * angle(cum_err);
%     
%     opt_sampling_inst = eps_hat(k) * conf.os_factor;
%     
%     switch interpolator_type
%         case 'none'
%             payload_data(k) = filtered_rx_signal(data_idx + round(opt_sampling_inst));
%             
%         case 'linear'
%             y = filtered_rx_signal(data_idx + floor(opt_sampling_inst) : data_idx + floor(opt_sampling_inst) + 1);
%             payload_data(k) = linear_interpolator(y, opt_sampling_inst - floor(opt_sampling_inst));
%             
%         case 'cubic'
%             y = filtered_rx_signal(data_idx + floor(opt_sampling_inst) - 1 : data_idx + floor(opt_sampling_inst) + 2);
%             payload_data(k) = cubic_interpolator(y, opt_sampling_inst - floor(opt_sampling_inst));
%             
%         otherwise
%             error('Unknown interpolator_type.');
%     end
%     

    payload_data(k) = filtered_rx_signal(data_idx);
    % Phase estimation    
    % Apply viterbi-viterbi algorithm
    deltaTheta = 1/4*angle(-payload_data(k)^4) + pi/2*(-1:4);
    
    % Unroll phase
    [~, ind] = min(abs(deltaTheta - theta_hat(k)));
    theta = deltaTheta(ind);
    
    % Lowpass filter phase
    theta_hat(k+1) = mod(0.01*theta + 0.99*theta_hat(k), 2*pi);
    
    % Phase correction
    payload_data(k) = payload_data(k) * exp(-1j * theta_hat(k+1));   % ...and rotate the current symbol accordingly
    
    data_idx = data_idx + conf.os_factor;
end

rxbits = payload_data;