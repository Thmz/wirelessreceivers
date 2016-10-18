clc % clear command window
clear % clear workspace

tic(); % Start time measurement
nbBits = 1000000; % number of bits 



%% Source: Generate random bits
txbits = randi([0 1],nbBits,1);



%% Mapping: Bits to symbols

% BPSK: if txbits = 1, symbol will be -1, if txbits = -1, symbol will be 1
tx = -(txbits *2  - 1); 

%% Channel: Apply BSC

flip = floor(rand(nbBits, 1) + 0.2); % generates vector with 80% zeros, 20% 1's
rx = tx .*( 1 - flip *2 ); % flip tx if flip = 1

%% Demapping: Symbols to bits

rxbits = (-rx + 1)/2;


%% BER: Count errors

cpr = rxbits ~= txbits;
errors = sum(cpr);


%% Output result
disp(['BER: ' num2str(errors/nbBits*100) '%'])

% Stop time measurement
toc()
