% Task 2.3.2
% Plot normalized correlator output 
%
%
plen        = 100;
SNR         = 55;

% generate praemble
preamble     = genpreamble(plen); 
bpskpreamble = -2*(preamble-0.5);

load ./task2.mat

% distort signal
rx_signal = signal + sqrt( 1 / 10^(SNR/10) /2) * (randn(size(signal))+1j*randn(size(signal)));


% correlate preamble with signal
c = detectpreamble(real(rx_signal),bpskpreamble,Inf);

% plot correlation sequence
figure(1)
plot(c)