% task_3_and_4 - Example solution for tasks 1.4.3 and 1.4.4
%    
% Author(s): Christian Senning, Nicholas Preyss, Alexios Balatsoukas-Stimming
%

clear all
close all
clc

% Initialization
SNR  = -6:2:12;
BERGray = zeros(size(SNR));
BERNonGray = zeros(size(SNR));

% number of bits
numSim = 10^4;

% Gray mapping
GrayMap = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];

% Non-Gray mapping
NonGrayMap = 1/sqrt(2) * [( 1-1j) ( 1+1j) (-1+1j) (-1-1j)];

% BER calculation for all SNR points
for ii = 1:numel(SNR)
   
    % Convert SNR from dB to linear
    SNRlin = 10^(SNR(ii)/10);
    
    % Generate source bitstream
    source = randi([0 1],numSim,2);
    
    % Map input bitstream using Gray mapping
    mappedGray = GrayMap(bi2de(source)+1).';
      
    % Add AWGN
    mappedGrayNoisy = mappedGray + sqrt(1/(2*SNRlin)) * (randn(size(mappedGray)) + 1i*randn(size(mappedGray)) );  
    
    % Demap
    [~,ind] = min((ones(numSim,4)*diag(GrayMap) - mappedGrayNoisy(:,[1 1 1 1])),[],2);
    demappedGray = de2bi(ind-1);
    
    % BER calculation for Gray mapping
    BERGray(ii) = mean(source(:) ~= demappedGray(:));
    
    % Map input bitstream using non-Gray mapping
    mappedNonGray = NonGrayMap(bi2de(source)+1).';
      
    % Add AWGN
    mappedNonGrayNoisy = mappedNonGray + sqrt(1/(2*SNRlin)) * (randn(size(mappedNonGray)) + 1i*randn(size(mappedNonGray)) );  
    
    % Demap
    [~,ind] = min((ones(numSim,4)*diag(NonGrayMap) - mappedNonGrayNoisy(:,[1 1 1 1])),[],2);
    demappedNonGray = de2bi(ind-1);
    
    % BER calculation for Gray mapping
    BERNonGray(ii) = mean(source(:) ~= demappedNonGray(:));

end

% graphical ouput
figure(3)
clf(3)
semilogy(SNR, BERGray, 'bx-' ,'LineWidth',3);
hold on
semilogy(SNR, BERNonGray, 'r*--','LineWidth',3);
xlabel('SNR (dB)')
ylabel('BER')
legend('Gray Mapping', 'Non-Gray Mapping')
grid on