clear
clc
% Oversampling factor
os_factor = 4;

% SNR
SNR = 3;
noAntenna = 4;

% transmitter
load task6_1


data_length = prod(image_size) * 8 / 2; % Number of QPSK data symbols
noframes = size(signal,1); 
symbolsperframe = data_length/noframes;

rxsymbols = zeros(noframes,symbolsperframe);

% Loop through all frames
for k=1:noframes
   
    Frame = signal(k,:);
    
    % Apply Rayleigh Fading Channel
    h = randn(noAntenna,1)+1i*randn(noAntenna,1);
    chanFrame = h * Frame;
    
    % Add White Noise
    SNRlin = 10^(SNR/10);
    noiseFrame = chanFrame + 1/sqrt(2*SNRlin)*(randn(size(chanFrame)) + 1i*randn(size(chanFrame)));

    
    
    % Matched Filter
    filtered_rx_signal = zeros(1);
    data_idx = 0;
    theta = 0;
    h = 0;
    
    for i = 1:noAntenna
        this_filtered_rx_signal = matched_filter(noiseFrame(i,:), os_factor, 6); % 6 is a good value for the one-sided RRC length (i.e. the filter has 13 taps in total) 
        % Frame synchronization
        [this_data_idx, this_theta, this_h] = frame_sync(this_filtered_rx_signal.', os_factor); % Index of the first data symbol
        
        % Selects signal with strongest channel
        if(abs(this_h) > abs(h))
            h = this_h;
            theta = this_theta;
            data_idx = this_data_idx;
            filtered_rx_signal = this_filtered_rx_signal/this_h;
        end
    end
    
    if(data_idx == 0)
        error('No syncronisation sequence found');
    end
    
    % Pick correct sampling points (no timing error)
    correct_samples = filtered_rx_signal(data_idx:os_factor:data_idx+(symbolsperframe*os_factor)-1);
    
    rxsymbols(k,:) = correct_samples;
    
end

combined_rxsymbols = reshape(rxsymbols.',1,noframes*symbolsperframe);

rxbitstream = demapper(combined_rxsymbols); % Demap Symbols
image_decoder(rxbitstream,image_size) % Decode Image