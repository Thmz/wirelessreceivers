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


% Downconvert
time = 1:1/conf.f_s:1+(length(rxsignal)-1)/conf.f_s;
rxsignal = rxsignal .* exp(-1j*2*pi*conf.f_c * time.');
rxsignal = lowpass(rxsignal, conf);

% Get datalength, depends on modulation type
switch(conf.modulation)
    case 'BSPK'
        data_length = conf.nbits;
    case 'QPSK'
        data_length = conf.nbits/2;
end;

% Receiver
filtered_rx_signal = matched_filter(rxsignal, conf.os_factor, 6); % 6 is a good value for the one-sided RRC length (i.e. the filter has 13 taps in total)

% Frame synchronization
[data_idx, theta] = frame_sync(filtered_rx_signal, conf.os_factor) % Index of the first data symbol
payload_data = zeros(data_length, 1); % The payload data symbols
theta_hat = zeros(data_length, 1);   % Estimate of the carrier phase
theta_hat(1) = theta;

% Loop over the data symbols with estimation and correction of phase
for k = 1 : data_length,

    % No time estimation needed due to very high oversampling factor.
    % Preamble detection will do sufficient time alignment
    
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
    
    % Move to next sample
    data_idx = data_idx + conf.os_factor;
end

% Demap
rxbits = demapper(payload_data,conf.modulation);