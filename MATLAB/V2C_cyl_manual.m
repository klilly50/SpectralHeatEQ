function vals = V2C_cyl_manual(vals,a)
% This code takes as argument a matrix of function values at n-dimensional
% Chebyshev and equidistant points, as well as an "order" string that tells
% the computer which dimensions are sampled at equidistant (Fourier)
% points. It cycles through each dimension, performing a DCT or DFT as
% appropriate.
% 
% In the "order" string, the character 't' (or 'T') denotes "trigonometric"
% dimensions. These are the ones for which the DFT is used. The character
% '0' is used to denote a dimension not to transform.
% 
% The time complexity of this algorithm is K*log(K), where K =
% prod(size(vals)). It is done in place, so space requirements are O(K).
% 
% August 21, 2017. David Darrow.

if ( size(a,1) > 1 )
    return
end

dim = size(a,2);
SIZE = size(vals);

if isstring(a)
    a = char(a);
end

for kk = 1:dim
    % First, "vectorize" to catch all columns at once.
    vals = reshape(vals,SIZE(kk),prod(SIZE)/SIZE(kk)); 
    
    % Transform values according to order.
    if ( a(kk) == 't' || a(kk) == 'T' )
        vals = trigv2c(vals);
    elseif ( a(kk) == '0' )
    else
        vals = chebv2c(vals);
    end
    
    if ( dim ~= 1 ) % It thinks 1D arrays are 2D, so "permute" causes problems
        vals = reshape(vals, [5 5 5]);
        vals = permute(vals, [2 3 1]);
    end

    if (a(kk) == 'r' || a(kk) == 'z')
        vals(:,:,2) = -vals(:,:,2);
        vals(:,:,4) = -vals(:,:,4);
    end

end

end