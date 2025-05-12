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
        vals = cheb2vec(vals); % run for equiv matlab code version
    end
    
    if ( dim ~= 1 ) % It thinks 1D arrays are 2D, so "permute" causes problems
        vals = reshape(vals,circshift(SIZE,dim-kk+1,2));
        vals = permute(vals,circshift(1:dim,dim-1,2));
    end
end

end

%% manual implementations
function coeffs = cheb2vec(vals)
    N = size(vals, 1);
    coeffs = zeros(size(vals));

    % Mirror and FFT each column
    for col = 1:size(vals,2)
        v = vals(:, col);
        v_ext = [v; v(end-1:-1:2)];
        c_ext = real(fft(v_ext));
        c = c_ext(1:N) / (N - 1);
        c(1) = c(1)/2;
        c(end) = c(end)/2;
        coeffs(:, col) = c;
    end
end

% incorrect
% function coeffs = trigv2c(vals)
%     N = length(vals);              % Number of sample points
%     coeffs = fftshift(fft(vals))/N; % Compute Fourier coefficients
% end

function coeffs = trigv2c(values)
%VALS2COEFFS   Convert values at N equally spaced points between [-1 1) 
%   to N trigonometric coefficients.
%   C = VALS2COEFFS(V) returns the vector of N coefficients such that: 
%   If N is odd
%       F(x) = C(1)*z^(-(N-1)/2) + C(2)*z^(-(N-1)/2-1) + ... + C(N)*z^((N-1)/2)
%   If N is even
%       F(x) = C(1)*z^(-N/2) + C(2)*z^(-N/2+1) + ... + C(N)*z^(N/2-1)           
%   where z = exp(1i*pi*x) and -1 <= x <= 1. 
%
%   F(x) interpolates the data [V(1) ; ... ; V(N)] at the N equally 
%   spaced points x_k = -1 + 2*k/N, k=0:N-1. 
%
%   If the input V is an (N+1)xM matrix, then C = VALS2COEFFS(V) returns the
%   (N+1)xM matrix of coefficients C such that F_j(x) intterpolates
%   [V(1,j) ; ... ; V(N+1,j)] for j=1:M using the same formula as above for
%   each column.
%
% See also COEFFS2VALS, TRIGPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% *Note about symmetry*.  Some of the code below is designed to
% enforce two symmetries whose failure might disturb users:
% VALUES exactly hermitian ==> COEFFS exactly real
% VALUES exactly skew-hermitian ==> % COEFFS exactly imaginary
% This is necessary because the MATLAB FFT code does not
% exactly preserve symmetries.

% Get the length of the input:
n = size(values, 1);

% Trivial case (constant or empty):
if ( n <= 1 )
    coeffs = values; 
    return
end

% test for symmetry
vals = double([values;values(1,:)]);
isHerm = max(abs(vals-conj(vals(end:-1:1, :))),[],1) == 0;
isSkew = max(abs(vals+conj(vals(end:-1:1, :))),[],1) == 0;

% compute coefficients
coeffs = (1/n)*fftshift(fft(values, [], 1), 1);

% correct if symmetric
coeffs(:,isHerm) = real(coeffs(:,isHerm));
coeffs(:,isSkew) = 1i*imag(coeffs(:,isSkew));

% These coefficients are for interpolation defined on [0,2*pi), but we want
% to work on [-pi,pi). To fix the coefficients for this we just need to
% assign c_k = (-1)^k c_k, for k=-(N-1)/2:(N-1)/2 for N odd and 
% k = -N/2:N/2-1 for N even.
if ( mod(n, 2) ) 
    even_odd_fix = (-1).^(-(n-1)/2:(n-1)/2).';
else
    even_odd_fix = (-1).^((-n/2):(n/2-1)).';
end

coeffs = bsxfun(@times, coeffs, even_odd_fix);

end
