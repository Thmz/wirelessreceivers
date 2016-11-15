clear
clc
% Oversampling factor
os_factor = 1;

% SNR
SNR = 50000;
noAntenna = 1;

% transmitter
load ./pn_sequence_fading.mat % loads as var signal, data symbols
load ./ber_pn_seq.mat % loads as ber_pn_seq, data bits

data_length = length(ber_pn_seq)/2; 
noframes = size(signal,1); % wil be 1
symbolsperframe = data_length/noframes;

rxsymbols = zeros(noframes,symbolsperframe);

% Loop through all frames
for k=1:noframes
   
    Frame = signal(k,:);
    
    % Apply Rayleigh Fading Channel
    h = randn(noAntenna,1)+1i*randn(noAntenna,1);
    chanFrame = Frame;
    
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
        [this_data_idx, this_theta, this_h] = frame_sync(this_filtered_rx_signal.', os_factor) % Index of the first data symbol
        
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

rxbits = demapper(combined_rxsymbols); % Demap Symbols

cpr = rxbits ~= ber_pn_seq;
errors = sum(cpr);
BER = errors/length(cpr);
disp(['BER: ' num2str(BER)])