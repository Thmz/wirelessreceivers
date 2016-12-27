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
tx = mapper(txbits, conf.modulation_order);
% length(tx_symbols) is in conf.data_length

%% Modify length of tx symbols
% Change into vector that is multiple of number of carriers

% The current size of the symbols is
conf.data_length;

% To be a multiple of the number of carriers, the required number of symbols is
tx_length = conf.n_carriers * conf.n_data_symbols;

% So we add more random symbols, to be exact amount of other random symbols is:
ext_length = tx_length - conf.data_length;

% Generate the symbols
ext_bits = randi([0 1],1, ext_length*2); % *2 because each symbol requires 2 bits
ext_symbols = mapper(ext_bits, 2);

% Add extension symbols to tx_symbols, so the length of tx_symbols is now a
% multiple of the number of carriers
tx = [tx(:); ext_symbols(:)];


%% Training frames
% A training symbol for each sub carrier frequency based on lfsr
% Training_symbols are the symbols for one training frame
[tx, nb_inserted] = insert_regular(tx, conf.training_symbols, conf.training_interval * conf.n_carriers);

% Update number of OFDM symbols because training symbols added
n_symbols = conf.n_data_symbols + nb_inserted; %  * length(conf.training_symbols) / conf.n_carriers;

%% Create a matrix of the symbols, with each column an OFDM symbol
% A bit inefficient but makes it easier to analyze frame data
txmatrix = reshape(tx, [conf.n_carriers, n_symbols]); % Matrix with OFDM symbols

%% Convert with IDFT operation
% and add cyclic prefix
symbol_length =  conf.n_carriers * conf.os_factor;
disp(['SYMBOL LENGTH ' num2str(symbol_length)])

s = zeros( (1+conf.cpref_length) * symbol_length, n_symbols);
for i = 1:n_symbols
   temp = osifft(txmatrix(:, i), conf.os_factor);
   cpref = temp(end - conf.cpref_length*symbol_length+1: end);
   s(:, i) = [cpref; temp];
end

%% Convert to serial
s = s(:);
disp(['LENGTH TX SIGNAL WITHOUT PREAMBLE ' num2str(length(s))])

%% Preamble
% Preamble: pulse shaping but no lowpass
preamble_bpsk = (1 - 2 * lfsr_framesync(conf.npreamble));
preamble_up = upsample(preamble_bpsk, conf.os_factor);
preamble_shaped = conv(preamble_up, conf.h.','same');

%% Normalize signal
figure, plot(abs([preamble_shaped ;s])), title('TX signal without normalization'), ylabel('amplitude'), xlabel('symbol nb');
normp = mean(abs(preamble_shaped).^2);
norms = mean(abs(s).^2);
txsignal = [ preamble_shaped /sqrt(normp); s/sqrt(norms)];
disp(['LENGTH TX SIGNAL WITH PREAMBLE ' num2str(length(txsignal))])
figure, plot(abs(txsignal)), title('TX signal with normalization'), ylabel('amplitude'), xlabel('symbol nb');


%% Apply lowpass
txsignal = ofdmlowpass(txsignal, conf, conf.corner_f);

%% Upconvert
time = 0:1/conf.f_s: (length(txsignal) -1)/conf.f_s;
txsignal = real(txsignal .* exp(1j*2*pi*conf.f_c * time.'));