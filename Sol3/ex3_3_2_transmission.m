function ex3_3_2_transmission
% ex3_3_2 - Task 3.3.2 solution
%    Measure BER of a root raised cosine shaped pulse transmission
%    
% 
% Author(s): Nicholas Preyss
% Copyright (c) 2012 TCL.
% 

SNR                 = 8;
rx_filterlen_array  = [0:10];
tx_filterlen        = 20; % tx_filterlen > rx_filterlen
os_factor           = 4; % oversampling factor
blank               = 10;

len          = 1000000;

% Generate random bitstream
bitstream = randi([0 1],1,len);

% Convert to QPSK symbol diracs
symbol1d            = 2*(bitstream-0.5);
symbol2s            = zeros(1,os_factor*(len/2));
symbol2s(1:4:end)   = 1/sqrt(2)*(symbol1d(1:2:end) + 1i*symbol1d(2:2:end)); 

% Create RRC pulse 
rolloff_factor = 0.22;
pulse = rrc(os_factor, rolloff_factor, tx_filterlen);

% Shape the symbol diracs with pulse
signal = conv([zeros(1,blank*os_factor) symbol2s zeros(1,blank*os_factor)],pulse.','same');

% add AWGN
rx_signal = awgn(signal,SNR);

% Simulate BER for each RX filter length
for i=1:length(rx_filterlen_array)
    rx_filterlen = rx_filterlen_array(i);
    
    filtered_rx_signal = matched_filter_np(rx_signal,os_factor,rx_filterlen);

    sampled_signal = filtered_rx_signal(1+blank*os_factor-rx_filterlen:4:end-blank*os_factor+rx_filterlen);

    decoded_bits = demapper(sampled_signal);

    ber(i) = sum(bitstream ~= decoded_bits.')/len;
end

% Plot results
figure(1)
semilogy(rx_filterlen_array,ber);
xlabel('Filter-length');
ylabel('SNR');
