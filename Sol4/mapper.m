function [txbits] = mapper(tx, modulation)

switch modulation
    case 'BPSK'
        txbits = 1-2*tx
    case 'QPSK'
        
       
end