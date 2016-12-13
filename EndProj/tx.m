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
tx_symbols = mapper(txbits, 2);
% length(tx_symbols) is in conf.data_length

% n_data_frames zit nu in config! onderstaande makes no sense 
% %% Convert to parallel
% n_frames = ceil(conf.data_length/conf.n_carriers)+1; % total number of OFDM symbols
% disp(['NUMBER OF SYMBOLS ' num2str(n_frames)])



%% Modify length of tx symbols
% Change into vector that is multiple of number of carriers

% The current size of the symbols is
conf.data_length;

% To be a multiple of the number of carriers, the required number of symbols is
tx_length = conf.n_carriers * conf.n_data_frames;

% So we add more random symbols, to be exact amount of other random symbols is:
ext_length = tx_length - conf.data_length

% Generate the symbols
ext_bits = randi([0 1],1, ext_length*2); % *2 because each symbol requires 2 bits
ext_symbols = mapper(ext_bits, 2);

%size1 = size(txsymbols(1:length(tx)+conf.n_carriers))
%size2 = size([training_symbols ;tx])

% Add extension symbols to tx_symbols, so the length of tx_symbols is now a
% multiple of the number of carriers
tx_symbols = [tx_symbols(:); ext_symbols(:)];


%% Training frames
% A training symbol for each sub carrier frequency based on lfsr
% Training_symbols are the symbols for one training frame
[tx_symbols, nb_inserted] = insert_regular(tx_symbols, conf.training_symbols, conf.training_interval * conf.n_carriers);

% Update number of frames because training symbols added
n_frames = conf.n_data_frames + nb_inserted; %  * length(conf.training_symbols) / conf.n_carriers;

%% Create a matrix of the symbols, with each column a frame
% A bit inefficient but makes it easier to analyze frame data
txmatrix = reshape(tx_symbols, [conf.n_carriers, n_frames]); % Matrix with OFDM symbols

%% Convert with IDFT operation
% and add cyclic prefix
symbol_length =  conf.n_carriers * conf.os_factor;
disp(['SYMBOL LENGTH ' num2str(symbol_length)])

s = zeros( (1+conf.cpref_length) * symbol_length, n_frames);
for i = 1:n_frames
   fttt = osifft(txmatrix(:, i), conf.os_factor);
   cpref = fttt(end - conf.cpref_length*symbol_length+1: end);
   s(:, i) = [cpref; fttt];
end

%% Convert to serial
s = s(:);

%% Calculate corner frequency
%corner_f = conf.f_bw*2; %???
%disp(['LPF CORNER FREQ ' num2str(corner_f)])
%s = ofdmlowpass(s, conf, corner_f);

disp(['LENGTH TX SIGNAL WITHOUT PREAMBLE ' num2str(length(s))])

%% Preamble
% Preamble: pulse shaping but no lowpass
preamble_bpsk = (1 - 2 * lfsr_framesync(conf.npreamble));
preamble_up = upsample(preamble_bpsk, conf.os_factor);
preamble_shaped = conv(preamble_up, conf.h.','same');

%% Get signal
normp = mean(abs(preamble_shaped).^2);
norms = mean(abs(s).^2);
txsignal = [ preamble_shaped /sqrt(normp); s/sqrt(norms)];
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