clear
clc

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
sigmaDeltaTheta = 0.001;
phase_noise = zeros(length(signal), 1);
phase_noise(1) = 2*pi*rand;
for i = 2:length(phase_noise)
    phase_noise(i) = mod(phase_noise(i-1) + sigmaDeltaTheta * randn(1, 1), 2*pi); 
end;

% Apply phase noise
antSig = antSig .* exp( 1i * phase_noise);

% Receiver
filtered_rx_signal = matched_filter(antSig, os_factor, 6); % 6 is a good value for the one-sided RRC length (i.e. the filter has 13 taps in total)

% Frame synchronization
[data_idx, theta] = frame_sync(filtered_rx_signal, os_factor); % Index of the first data symbol
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
   theta = 1/4 * angle( - payload_data(k)^4)  +  pi/2 * (1:4);
    
    % Phase correction
   aa = exp(1i*theta);
   bb = exp(1i*theta_hat(k));
   dist = (real(aa) - real(bb)).^2 + (imag(aa) - imag(bb)).^2;
    
   [M, I] = min(dist);
   theta_closest = theta(I);
   theta_hat(k+1) = mod(theta_closest, 2*pi);
    
   payload_data(k) = payload_data(k) * exp(-1i * theta_hat(k+1));

    
    data_idx = data_idx + os_factor;
end

% Plot phase
figure(1),
plot(theta_hat)
hold on
plot(phase_noise(data_idxInit:os_factor:end), '-r')
a = axis;
a(3:4) = [0,2*pi];
axis(a)
set(gca,'ytick',0:pi/4:2*pi)
grid on;

% Draw image
image_decoder(demapper(payload_data), image_size);