function [txsignal conf] = tx(txbits,conf,k)
% Digital Transmitter
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete transmitter
%   consisting of:
%       - modulator
%       - pulse shaping filter
%       - up converter
%   in digital domain.
%
%   txbits  : Information bits
%   conf    : Universal configuration structure
%   k       : Frame index
%


% dummy 400Hz sinus generation
%time = 1:1/conf.f_s:4;
%txsignal = 0.3*sin(2*pi*400 * time.');

preamble_bpsk = (1 - 2 * lfsr_framesync(conf.npreamble));
tx = [preamble_bpsk; mapper(txbits, conf.modulation_order)];

%oversample
tx = upsample(tx, conf.os_factor);

% Shape the symbol diracs with pulse
txsignal = conv(tx, conf.h.','same');

% Upconvert
time = 0:1/conf.f_s: (length(txsignal) -1)/conf.f_s;
txsignal = real(txsignal .* exp(1j*2*pi*conf.f_c * time.'));