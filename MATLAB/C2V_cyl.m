function vals = C2V_cyl(vals,a)
% This code takes as argument a matrix of coefficients for an
% n-dimensional Chebyshev/Fourier series, as well as an "order" string that
% tells the computer which dimensions use Fourier coefficients. It cycles 
% through each dimension, performing an iDCT or iDFT as appropriate.
% 
% In the "order" string, the character 't' (or 'T') denotes "trigonometric"
% dimensions. These are the ones for which the iDFT is used. The character
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
    vals = reshape(vals,SIZE(kk),prod(SIZE)/SIZE(kk));
    if ( a(kk) == 't' || a(kk) == 'T' )
        for ii = 1:prod(SIZE)/SIZE(kk)
            vals(:,ii) = trigtech.coeffs2vals(vals(:,ii));
        end
    elseif ( a(kk) == '0' )
    else
        for ii = 1:prod(SIZE)/SIZE(kk)
            vals(:,ii) = chebtech2.coeffs2vals(vals(:,ii));
        end
    end

    if ( dim ~= 1 ) % It thinks 1D arrays are 2D, so "permute" causes problems
        vals = reshape(vals,circshift(SIZE,dim-kk+1,2));
        vals = permute(vals,circshift(1:dim,dim-1,2));
    end

    fprintf('dimension %s\n', a(kk));
    disp(vals);

end

vals = real(vals);
end