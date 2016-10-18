clear 
clc
close all

profile on

%% Generate preamble
Np = 100; % length of preamble
treshold = 15; % Treshold for peak detection, 15 is good value, 3 is bad value


reg = ones(1, 8);
pbits = zeros(1,Np); % output

for i = 1:Np
    
    pbits(i) = reg(8);
    reg(2:8) = reg(1:7);
    reg(1) = xor(xor(xor(pbits(i), reg(7)),reg(6)),reg(5)) ;
    
end
pbits;

%% Mapping: Bits to symbols

% Preamble in BPSK: if txbits = 1, symbol will be -1, if txbits = -1, symbol will be 1
p = -(pbits * 2  - 1); 


%% Received signal
%Load mapped signal
load task2.mat
tx = signal;

% According to assignment no noise
% SNR_dB = 10;
% wx = gen_noise(SNR_dB, tx);
% rx = tx + wx;
rx = tx;

correlated = zeros(length(rx), 1);  

%% Loop over all incoming symbols to calculate correlation and normalize
for n = Np:length(rx) 
    rx_part = rx(n-Np+1:n);

    % Correlate
    correlated(n) = conj(p) * rx_part;

    % Normalize
    normalized = abs(correlated(n)).^2/sum(abs(rx_part).^2);
    if(normalized > treshold)
        % Get image signal and display image
        image = rx(n+1:n+ (image_size(1) * image_size(2)) *4);
        imagebits = demapper(image);
        image_decoder(imagebits, image_size);
        break;
    end
end

profile off
profsave
profile viewer