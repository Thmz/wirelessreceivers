function [ out, n_parts ] = insert_regular( vec, insert_symbols, insert_interval )
%INSERT_REGULAR Inserts the given vector sequence insert_symbols at regular intervals in
%the given vector vec. The number of symbols between the inserted symbols will be the
%integer insert_interval
% Output is the vector vec with the inserted symbols.

% Make the length of the given vector a multiple of insert_interval
n_parts = ceil(length(vec) / insert_interval); % Will be the number of inserted training frames
veczeros = zeros(1, n_parts * insert_interval);
veczeros(1:length(vec)) = vec;

% Make a matrix of vec
vecmatrix = reshape(veczeros, [insert_interval, n_parts]);

% Make a matrix of training symbols
insert_symbols = insert_symbols(:); % Make sure it is a column vector
insmatrix = repmat(insert_symbols, 1, n_parts);

% Add insert symbols
outmatrix = [insmatrix; vecmatrix] ;

% Convert to column vector and remove unnecessary zero symbols at the end
out = outmatrix(1: (length(vec) + length(insert_symbols)*n_parts)); 

end

