function [ c ] = correlator( p, rx)
    c = zeros(length(rx), 1);

    for n = length(p):length(rx);
        c(n) = conj(p) * rx(n-length(p)+1:n);
    end
%     c = xcorr(p, rx);
%     c = c(length(rx):-1:1);
end

