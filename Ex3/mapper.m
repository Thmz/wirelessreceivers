function [ qpsk ] = mapper( bits )
% MAP map bit stream to QPSK
%   
qpsk = ((bits(1:2:length(bits)) - 0.5) + (bits(2:2:length(bits)) - 0.5)*1i ) * sqrt(2);

end

