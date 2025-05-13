% manual implementation of chebtech2.vals2coeffs

function coeffs = chebv2c(vals)
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