% Oversampling factor
os_factor = 4;

% SNR
SNR = 12;

% Interpolator type
interpolator_type = 'linear';

% transmitter
load task4
data_length = prod(image_size) * 8 / 2; % Number of QPSK data symbols

% Channel
antSig = awgn(signal, SNR);

% Create phase noise
theta_n = zeros(length(antSig),1);
theta_n(1) = 2*pi*rand(1);
sigmaDeltaTheta = 0.001;
for i = 2 : length(antSig)
   theta_n(i) = mod(theta_n(i-1) + sigmaDeltaTheta*randn,2*pi);
end

% Apply phase noise
antSig = antSig.*exp(1j*theta_n);

% Receiver
filtered_rx_signal = matched_filter(antSig, os_factor, 6); % 6 is a good value for the one-sided RRC length (i.e. the filter has 13 taps in total)

% Frame synchronization
[data_idx theta] = frame_sync(filtered_rx_signal, os_factor); % Index of the first data symbol
data_idxInit = data_idx;

payload_data = zeros(data_length, 1); % The payload data symbols
eps_hat = zeros(data_length, 1);     % Estimate of the timing error
theta_hat = zeros(data_length, 1);   % Estimate of the carrier phase
theta_hat(1) = theta;

% Loop over the data symbols with estimation and correction of the timing
% error and phase
for k = 1 : data_length,
    
    % timing estimator
    eps_hat(k) = timing_estimator(filtered_rx_signal(data_idx : data_idx + os_factor - 1));
    opt_sampling_inst = eps_hat(k) * os_factor;
    
    switch interpolator_type
        case 'none'
            payload_data(k) = filtered_rx_signal(data_idx + round(opt_sampling_inst));
            
        case 'linear'
            y = filtered_rx_signal(data_idx + floor(opt_sampling_inst) : data_idx + floor(opt_sampling_inst) + 1);
            payload_data(k) = linear_interpolator(y, opt_sampling_inst - floor(opt_sampling_inst));
            
        case 'cubic'
            y = filtered_rx_signal(data_idx + floor(opt_sampling_inst) - 1 : data_idx + floor(opt_sampling_inst) + 2);
            payload_data(k) = cubic_interpolator(y, opt_sampling_inst - floor(opt_sampling_inst));
            
        otherwise
            error('Unknown interpolator_type.');
    end
    
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
    
    data_idx = data_idx + os_factor;
end

% Plot phase
figure(1),
plot(theta_hat)
hold on
plot(theta_n(data_idxInit:os_factor:end), '-r')
a = axis;
a(3:4) = [0,2*pi];
axis(a)
set(gca,'ytick',0:pi/4:2*pi)
grid on;

% Draw image
image_decoder(demapper(payload_data), image_size);