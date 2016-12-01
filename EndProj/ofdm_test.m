
 % Generate random data
txbits = randi([0 1],conf.nbits,1);

%% Convert to QPSK
txbits = mapper(txbits, conf.modulation_order)];

%% Convert to parallel

n_samples = ceil(length(txbits)/conf.n_carriers)

length(txbits)
%a = txbits(1:256)
txzeros = zeros(1, conf.n_carriers*n_samples);
txzeros(1:length(txbits)) = txbits;
length(txzeros)

a = zeros(conf.n_carriers, n_samples);

for i = 1:n_samples
    a(:, i) = txzeros((i-1)*conf.n_carriers+1: i*conf.n_carriers);
end

%% Convert with IDFT operation

for i:n_samples
    s(i) = OSIFFT(a(i), conf.os_factor);
end

%% Convert to serial


% 

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