function ex2_3a_detectframe()
% ex2_4 - Example solution for task A
%    Display Image
%    
% 
% Author(s): Christian Senning, Nicholas Preyss
% Copyright (c) 2012 TCL.

plen        = 100;
SNR         = 10;
threshold   = 50;

% generate praemble
preamble     = genpreamble(plen); 
bpskpreamble = -2*(preamble-0.5);

load ./task2.mat
imagelen = prod(image_size)*4;

% distort signal
rx_signal = signal + sqrt( 1 / 10^(SNR/10) /2) * (randn(size(signal))+1j*randn(size(signal)));

% correlate preamble with signal
start = detectpreamble(real(rx_signal),bpskpreamble,threshold);

image_decoder(demapper(rx_signal(start:start+imagelen-1)), image_size);

end