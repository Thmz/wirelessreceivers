function image = bi2de(b)

% Converts a binary row vector b to a nonnegative decimal integer. 

% If b is a matrix, each row is interpreted separately as a binary number. 
% In this case, the output d is a column vector, each element of which is 
% the decimal representation of the corresponding row of b.

% Note that bi2de interprets the first column of b as the lowest-order digit.

image = b * 2.^(0:size(b,2)-1)';

end