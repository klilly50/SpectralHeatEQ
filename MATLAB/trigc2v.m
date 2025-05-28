function values = trigc2v(coeffs)
    %COEFFS2VALS   Convert trigonometric coefficients to values on [-1,1)
    %   values = trigc2v(coeffs) computes function values at equispaced
    %   points x_k = -1 + 2*k/N for k = 0:N-1, from the given trigonometric
    %   coefficients.
    %
    %   This is the inverse of trigv2c.
    
    n = size(coeffs, 1);
    
    % Fix sign alternation: invert the earlier transformation
    if mod(n, 2)
        even_odd_fix = (-1).^(-(n-1)/2:(n-1)/2).';
    else
        even_odd_fix = (-1).^(-n/2:n/2 - 1).';
    end
    coeffs = bsxfun(@times, coeffs, even_odd_fix);
    
    % Invert FFT shift and scale
    values = ifft(ifftshift(coeffs, 1), [], 1) * n;
    
    % Ensure Hermitian symmetry is respected for real outputs
    % Force real or imaginary as needed (optional, for robustness)
    % vals_full = [values; values(1,:)];
    % isHerm = max(abs(vals_full - conj(vals_full(end:-1:1,:))), [], 1) == 0;
    % isSkew = max(abs(vals_full + conj(vals_full(end:-1:1,:))), [], 1) == 0;
    % 
    % values(:, isHerm) = real(values(:, isHerm));
    % values(:, isSkew) = imag(values(:, isSkew)) * 1i;
end

