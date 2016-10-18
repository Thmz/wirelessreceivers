clc % clear command window
clear % clear workspace

nbBits = 1000000; % number of bits 
os_factor = 4;
tx_filtertaps = 41;
rx_filterlen = 1;
filter_alpha = 0.22;
SNR_dB = 8;
SNR = 10^(SNR_dB/10);

%% Source: Generate random bits
txbits = randi([0 1],nbBits,1);

%% Map to QPSK
tx = mapper(txbits);
tx_up = upsample(tx, os_factor);
% or, without using the upsample function:
% tx_up = zeros(length(tx)*4, 1);
% tx_up(1:os_factor:length(tx)*4) = tx ;

%% Shape RRC pulses
fir_rrc =  rrc(os_factor, filter_alpha, (tx_filtertaps - 1)/2 );
tx_up_wave = conv(tx_up, fir_rrc, 'same');

%% Add AWGN noise
wx =  1/sqrt(SNR*2) * (randn(length(tx_up_wave), 1) + randn(length(tx_up_wave), 1)*1i);
rx_wave = wx + tx_up_wave;

for rx_filterlen = 1:15
    
    %% Apply MF
    fir_rrc2 =  rrc(os_factor, filter_alpha, rx_filterlen);
    rx = conv(rx_wave, fir_rrc2, 'same');

    %% Downsample
    rx_down = downsample(rx, os_factor);
    % or without downsample function
    % rx_down = rx(1:os_factor:end);

    %% Demap
    rxbits = demapper(rx_down);

    %% Calculate BER
    BER = calc_ber(txbits, rxbits);
    disp(['FILTER TAPS: ' num2str(rx_filterlen*2+1) ', FILTER LENGTH: ' num2str(rx_filterlen)  ', BER: ' num2str(BER) ''])
end;
