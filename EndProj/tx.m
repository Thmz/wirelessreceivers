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

%% Convert to QPSK
tx = mapper(txbits, 2);
% length(tx) is in conf.data_length

%% Convert to parallel
n_symbols = ceil(conf.data_length/conf.n_carriers); % total number of OFDM symbols
disp(['NUMBER OF SYMBOLS ' num2str(n_symbols)])

% Change into vector that is multiple of number of carriers
txzeros = zeros(1, conf.n_carriers*n_symbols);
txzeros(1:length(tx)) = tx;
txmatrix = reshape(txzeros, [conf.n_carriers, n_symbols]); % Matrix with OFDM symbols

%% Convert with IDFT operation
% and add cyclic prefix
symbol_length =  conf.n_carriers * conf.os_factor;
disp(['SYMBOL LENGTH ' num2str(symbol_length)])

s = zeros( (1+conf.cpref_length) * symbol_length, n_symbols);
for i = 1:n_symbols
   fttt = osifft(txmatrix(:, i), conf.os_factor);
   cpref = fttt(end - conf.cpref_length*symbol_length+1: end);
   s(:, i) = [cpref; fttt];
end

%% Convert to serial
s = s(:)*1000;

%% Apply low pass
corner_f = conf.f_bw*2;
disp(['LPF CORNER FREQ ' num2str(corner_f)])
%s = ofdmlowpass(s, conf, corner_f);

disp(['LENGTH TX SIGNAL WITHOUT PREAMBLE ' num2str(length(s))])

%% Preamble
% Preamble: pulse shaping but no lowpass
preamble_bpsk = (1 - 2 * lfsr_framesync(conf.npreamble));
preamble_up = upsample(preamble_bpsk, conf.os_factor);
preamble_shaped = conv(preamble_up, conf.h.','same');

%% Get signal
txsignal = [ preamble_shaped; s];

disp(['LENGTH TX SIGNAL WITH PREAMBLE ' num2str(length(txsignal))])

%% Upconvert
time = 0:1/conf.f_s: (length(txsignal) -1)/conf.f_s;
txsignal = real(txsignal .* exp(1j*2*pi*conf.f_c * time.'));


%% TODO upconvert to carrier freq

% Old TX
% 
% 
% 
% for i = 1:conf.n_carriers:length(txbits)-conf.n_carriers
%     a(i, :) = txbits(i:i+conf.n_carriers-1);
% end
% a
%     
% 
% 
% 
% % 
% 
% 
% 
% 
% 
% 
% 
% % dummy 400Hz sinus generation
% %time = 1:1/conf.f_s:4;
% %txsignal = 0.3*sin(2*pi*400 * time.');
% 
% preamble_bpsk = (1 - 2 * lfsr_framesync(conf.npreamble));
% tx = [preamble_bpsk; mapper(txbits, conf.modulation_order)];
% 
% %oversample
% tx = upsample(tx, conf.os_factor);
% 
% % Shape the symbol diracs with pulse
% txsignal = conv(tx, conf.h.','same');
% 
% % Upconvert
% time = 0:1/conf.f_s: (length(txsignal) -1)/conf.f_s;
% txsignal = real(txsignal .* exp(1j*2*pi*conf.f_c * time.'));