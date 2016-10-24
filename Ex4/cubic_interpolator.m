function [ out ] = cubicinterpolator( f, k, delta )
    assert (delta > 0)
    
    A = [-1  3 -3  1 
          3 -6  3  0
         -2 -3  6 -1
          0  6  0  0];
    ff = [ f(k-1) f(k) f(k+1) f(k+2) ].';
    
    
    abcd = 1/6 * A * ff;
    
    x = delta;
    xx = [ x^3 x^2 x 1 ].';
    
    out = sum(abcd .* xx);
end


