function [rxbits] = demapper(rx, modulation)

switch modulation
    case 1 % BPSK
        rxbits = 1-(rx > 0);
    case 2 % QPSK
        
        % Convert noisy QPSK symbols into a bit vector. Hard decisions.
        
        a = rx(:); % Make sure "a" is a column vector
        
        b = [real(a) imag(a)] > 0;
        
        % Convert the matrix "b" to a vector, reading the elements of "b" rowwise.
        b = b.';
        b = b(:);
        
        rxbits = double(b); % Convert data type from logical to double
end