clear
clc
close all

%% Generate preamble
Np = 100; % length of preamble

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
tx = signal; %tx = signal(8000:12000); % to limit signal length



SNR_values = [-5 0 5 10]
for SNR_dB = SNR_values

    wx = gen_noise(SNR_dB, tx);
    rx = tx + wx;

    %% Calculate correlation
    correlated = correlator(p, rx);

    %% Normalize
    normalized = zeros(length(rx), 1);

    a = abs(correlated).^2;
    for n = Np:length(correlated)
        normalized(n) = a(n)/sum(abs(rx(n-Np+1:n)).^2);
    end


    %% Plot
    treshold = 15;

    figure, plot(abs(correlated)) % downsample before plotting, due to memory usage
    xlabel('symbol');
    ylabel('correlation');
    title(strcat('Correlation for Noise SNR =', num2str(SNR_dB), ' dB'));
    figure, plot(normalized);
    hold on;
    plot(treshold*ones(length(normalized), 1)); % downsample before plotting, due to memory usage
    xlabel('symbol');
    ylabel('normalized value');
    title(strcat('Normalized correlator output for Noise SNR =', num2str(SNR_dB), ' dB'));

end