N = 1000000;
p = 0.2;

% Start time measurement
tic();

% Source: Generate random bits
txbits = randi([0 1],N,1) ; 

% Mapping: Bits to symbols
tx = 2*txbits - 1;

% Channel: Apply binary symmetric channel
flip = -1*( 2*(rand(N,1) < p) - 1 ) ;
rx = tx .* flip;

% Demapping: Symbols to bits
rxbits = (rx == 1);

% BER: Count errors
errors = sum(rxbits ~= txbits);

% Stop time measurement
toc()

% Output result
disp(['BER: ' num2str(errors/length(rxbits)*100) '%'])
