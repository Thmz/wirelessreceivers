function [ out ] = linear_interpolator( f, k, delta )
    assert (delta > 0)
    out = f(k)*(1-delta) + f(k+1)*delta;
end

