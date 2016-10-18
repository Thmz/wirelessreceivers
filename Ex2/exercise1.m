Np = 100 % length of preamble

reg = ones(1, 8);
p = zeros(1,Np); % output

for i = 1:Np
    
    p(i) = reg(8);
    reg(2:8) = reg(1:7);
    reg(1) = mod(reg(8) + reg(4)+ reg(5) + reg(6),2) ;
    
end

p