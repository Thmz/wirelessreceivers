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

disp('****** START RX ******')



% Downconvert
time = 0:1/conf.f_s: (length(rxsignal) -1)/conf.f_s;
rxsignal = rxsignal .* exp(-1j*2*pi*(conf.f_c + conf.offset) * time.');
corner_f = conf.f_bw*1.1; % TODO *2 is te hoog wss
disp(['LPF CORNER FREQ ' num2str(corner_f)]);
rxsignal = ofdmlowpass(rxsignal, conf, corner_f);



%% Frame sync
[data_idx, ~] = frame_sync(rxsignal, conf.os_factor); % Index of the first data symbol
disp(['DATA IDX ' num2str(data_idx)]);


os_frame_length = conf.n_carriers*conf.os_factor;
n_data_frames = ceil(conf.data_length/conf.n_carriers);
n_training_frames = ceil(n_data_frames/conf.training_interval);
n_frames = n_data_frames + n_training_frames;

signal_length = (1+conf.cpref_length)*os_frame_length*n_frames; % length of total signal

disp(['FRAME LENGTH WITHOUT PREFIX ' num2str(os_frame_length)]);
disp(['NUMBER OF FRAMES ' num2str(n_frames)]);
disp(['EXPECTED SIGNAL LENGTH ' num2str(signal_length)]);

rxsignal = rxsignal(data_idx:data_idx+signal_length-1); % The payload data symbols




%% Serial to parallel
rxmatrix = reshape(rxsignal, [ (1+conf.cpref_length)*os_frame_length  n_frames]);

%% Delete cyclic prefix
rxmatrix = rxmatrix(conf.cpref_length*conf.n_carriers*conf.os_factor+1: end, :);

%% FFT
rx = zeros(conf.n_carriers, n_frames);
for i = 1:n_frames
    rx(:, i) = osfft(rxmatrix(:,i), conf.os_factor);
end

%% Process signal
% Can also be done in loop above

payload_data = [];
theta_hat = zeros(conf.n_carriers, n_frames+1);

% Loop over frames
for i = 1:n_frames
    frame_data = rx(:, i);
    
    % Check if training frame
    is_training = mod(i-1, conf.training_interval+1) == 0;
    
    if(is_training)
        %% Check phase of training bits
        training_bits_phase = mod(angle(frame_data),2*pi);
        orig_training_bits_phase = mod(angle((1 - 2 * lfsr_framesync(conf.n_carriers))), 2*pi);
        phase_difference =  training_bits_phase - orig_training_bits_phase ;
        theta_hat(:, i+1) = phase_difference;
    else
        
        if(conf.do_phase_estim)
            if(conf.do_phase_track)
                
                % Apply Viterbi-Viterbi phase estimation
                deltaTheta = 1/4*angle( repmat(-frame_data(:).^4, 1 , 6)) + repmat( pi/2*(-1:4), conf.n_carriers, 1);
                
                % Unroll phase
          
                [~, ind] = min(abs(deltaTheta - repmat(theta_hat(:,i),1, 6)), [], 2);
                indvec = (0:conf.n_carriers-1).*6 + ind'; 
                deltaTheta = deltaTheta';
                theta = deltaTheta(indvec)';
                
                % Lowpass filter phase
                theta_hat(:, i+1) = mod(0.01*theta + 0.99*theta_hat(:, i), 2*pi);
            else
                theta_hat(:,i+1) = theta_hat(:, i);
            end
            % Phase correction
            frame_data = frame_data .* exp(-1j * theta_hat(:, i+1));   % ...and rotate the current symbol accordingly
        end        
        % Add values to received payload data
        payload_data = [payload_data(:) ; frame_data(:)];
    end  
end
figure, plot(theta_hat(1:10, :)')
title('Theta hat');


%phase_difference;

% %% Correct phase
% for i = 0:n_frames-2
%     sizerx = size(rx(i*conf.n_carriers +1 : conf.n_carriers*(i+1)))
%    rx(i*conf.n_carriers +1 : conf.n_carriers*(i+1)) = rx(i*conf.n_carriers +1 : conf.n_carriers*(i+1)) .*exp(-j*phase_difference);
% end

%% Demap
rxbits = demapper(payload_data, 2);

%% Delete unnecessary data at the end
rxbits = rxbits(1:conf.data_length*2);






%% OLD

% % Downconvert
% time = 1:1/conf.f_s:1+(length(rxsignal)-1)/conf.f_s;
% rxsignal = rxsignal .* exp(-1j*2*pi*(conf.f_c + conf.offset) * time.');
% rxsignal = lowpass(rxsignal, conf);
%
% % Receiver pulse shaping
% filtered_rx_signal = conv(rxsignal, conf.h.','same');
%
% % Frame synchronization
% [data_idx, theta] = frame_sync(filtered_rx_signal, conf.os_factor) % Index of the first data symbol
% payload_data = zeros(conf.data_length, 1); % The payload data symbols
% theta_hat = zeros(conf.data_length, 1);   % Estimate of the carrier phase
% theta_hat(1) = theta;
%
% % Loop over the data symbols with estimation and correction of phase
% for k = 1 : conf.data_length,
%
%     % No time estimation needed due to very high oversampling factor.
%     % Preamble detection will do sufficient time alignment
%
%     payload_data(k) = filtered_rx_signal(data_idx);
%
%     % Phase estimation
%     % Apply viterbi-viterbi algorithm
%     deltaTheta = 1/4*angle(-payload_data(k)^4) + pi/2*(-1:4);
%
%     % Unroll phase
%     [~, ind] = min(abs(deltaTheta - theta_hat(k)));
%     theta = deltaTheta(ind);
%
%     % Lowpass filter phase
%     theta_hat(k+1) = mod(0.01*theta + 0.99*theta_hat(k), 2*pi);
%
%     % Phase correction
%     payload_data(k) = payload_data(k) * exp(-1j * theta_hat(k+1));   % ...and rotate the current symbol accordingly
%
%     % Move to next sample
%     data_idx = data_idx + conf.os_factor;
% end
%
% % Demap
% rxbits = demapper(payload_data,conf.modulation_order);