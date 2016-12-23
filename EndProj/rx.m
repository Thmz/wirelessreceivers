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
plot(abs(fft(rxsignal)))
time = 0:1/conf.f_s: (length(rxsignal) -1)/conf.f_s;
rxsignal = rxsignal .* exp(-1j*2*pi*(conf.f_c + conf.offset) * time.');
rxsignal = ofdmlowpass(rxsignal, conf, conf.corner_f);



%% Frame sync
[data_idx, ~] = frame_sync(rxsignal, conf.os_factor); % Index of the first data symbol
disp(['DATA IDX ' num2str(data_idx)]);


os_symbol_length = conf.n_carriers*conf.os_factor;
n_data_symbols = ceil(conf.data_length/conf.n_carriers);
n_training_symbols = ceil(n_data_symbols/conf.training_interval);
n_symbols = n_data_symbols + n_training_symbols;

frame_length = (1+conf.cpref_length)*os_symbol_length*n_symbols; % length of total signal

disp(['SYMBOL LENGTH WITHOUT PREFIX ' num2str(os_symbol_length)]);
disp(['NUMBER OF SYMBOLS ' num2str(n_symbols)]);
disp(['EXPECTED FRAME LENGTH ' num2str(frame_length)]);

rxsignal = rxsignal(data_idx:data_idx+frame_length-1); % The payload data symbols




%% Serial to parallel
rxmatrix = reshape(rxsignal, [ (1+conf.cpref_length)*os_symbol_length  n_symbols]);

%% Delete cyclic prefix
rxmatrix = rxmatrix(conf.cpref_length*conf.n_carriers*conf.os_factor+1: end, :);

%% FFT
rx = zeros(conf.n_carriers, n_symbols);
for i = 1:n_symbols
    rx(:, i) = osfft(rxmatrix(:,i), conf.os_factor);
end

%% Process signal
% Can also be done in loop above

payload_data = [];
theta_hat = zeros(conf.n_carriers, n_symbols+1); % without phase tracking
theta_hat_pt = zeros(conf.n_carriers, n_symbols+1); % With phase tracking, doing both to compare

% Loop over frames
for i = 1:n_symbols
    symbol_data = rx(:, i);
    
    % Check if training frame
    is_training = mod(i-1, conf.training_interval+1) == 0;
    
    if(is_training)
        channel = symbol_data./conf.training_symbols;
        figure('Name', 'channel magnitude'), plot(abs(channel));
        figure('Name', 'channel phase'), plot(angle(channel));
        figure('Name', 'channel impulse response'), plot(abs(ifft(channel)));
        
        % Check phase of training bits
        training_bits_phase = mod(angle(symbol_data),2*pi);
        orig_training_bits_phase = mod(angle((1 - 2 * lfsr_framesync(conf.n_carriers))), 2*pi);
        phase_difference =  training_bits_phase - orig_training_bits_phase ;
        theta_hat(:, i+1) = phase_difference;
        theta_hat_pt(:, i+1) = phase_difference;
    else
       % symbol_data = symbol_data ./ channel; % Rescale symbols
        
        if(conf.do_phase_estim)
          
                
            % Apply Viterbi-Viterbi phase estimation
            deltaTheta = 1/4*angle( repmat(-symbol_data(:).^4, 1 , 6)) + repmat( pi/2*(-1:4), conf.n_carriers, 1);

            % Unroll phase
            [~, ind] = min(abs(deltaTheta - repmat(theta_hat_pt(:,i),1, 6)), [], 2);
            indvec = (0:conf.n_carriers-1).*6 + ind'; 
            deltaTheta = deltaTheta';
            theta = deltaTheta(indvec)';
            % Lowpass filter phase
            theta_hat_pt(:, i+1) = mod(0.05*theta + 0.95*theta_hat_pt(:, i), 2*pi);
             
            % For the case without phase tracking
            theta_hat(:, i+1) = theta_hat(:, i);
            if(conf.do_phase_track)    
               theta_corr = theta_hat_pt(:, i+1); 
            else
               theta_corr = theta_hat(:, i+1);
            end
            % Phase correction
            symbol_data = symbol_data .* exp(-1j * theta_corr);   % ...and rotate the current symbol accordingly
        end        
        % Add values to received payload data
        payload_data = [payload_data(:) ; symbol_data(:)];
    end  
end

% Plot 30 carriers and their estimation of theta
figure, plot(theta_hat(1:20:256, :)')
title('Estimation of theta (Theta hat) without phase tracking');


figure, plot(theta_hat_pt(1:20:256, :)')
title('Estimation of theta (Theta hat) with phase tracking');


%% Demap
rxbits = demapper(payload_data, conf.modulation_order);

%% Delete unnecessary data at the end
rxbits = rxbits(1:conf.data_length*2);
