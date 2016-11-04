% Oversampling factor
os_factor = 4;

% SNR
SNR = 6;
noAntenna = 3;

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
    plot(chanFrame)
    %
    % Receiver with Single Antenna
    %
    
    % Matched Filter
    
    filtered_rx_signal = matched_filter(noiseFrame(1,:), os_factor, 6); % 6 is a good value for the one-sided RRC length (i.e. the filter has 13 taps in total)

    % Frame synchronization
    [data_idx theta magn] = frame_sync(filtered_rx_signal.', os_factor); % Index of the first data symbol
    
    % Pick correct sampling points (no timing error)
    correct_samples = filtered_rx_signal(data_idx:os_factor:data_idx+(symbolsperframe*os_factor)-1);
    
    rxsymbols(k,:) = correct_samples;
    
end

combined_rxsymbols = reshape(rxsymbols.',1,noframes*symbolsperframe);

rxbitstream = demapper(combined_rxsymbols); % Demap Symbols
image_decoder(rxbitstream,image_size) % Decode Image