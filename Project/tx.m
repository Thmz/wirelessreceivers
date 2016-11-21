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

framesync_bpsk = (1 - 2 * lfsr_framesync(conf.npreamble));
tx = [zeros(100, 1); framesync_bpsk; mapper(txbits, conf.modulation)];

%oversample
tx = upsample(tx, conf.os_factor);

% Create RRC pulse
tx_filterlen = 20;
rolloff_factor = 0.22;
pulse = rrc(conf.os_factor, rolloff_factor, tx_filterlen);

% Shape the symbol diracs with pulse
txsignal = conv(tx, pulse.','same');

% Upconvert
time = 1:1/conf.f_s:1+ (length(txsignal) -1)/conf.f_s;
txsignal = real(txsignal .* exp(1j*2*pi*conf.f_c * time.'));
