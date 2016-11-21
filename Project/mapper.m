function [tx] = mapper(txbits, modulation)

switch modulation
    case 'BPSK'
        tx = 1-  2 * txbits;
    case 'QPSK'
       symbol1d            = 2*(txbits-0.5);
      % symbol2s            = zeros(1,os_factor*(len/2));
       tx   = 1/sqrt(2)*(symbol1d(1:2:end) + 1i*symbol1d(2:2:end)); 
end