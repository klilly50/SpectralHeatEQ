%% testing FFT function in MATLAB

f = @(x) exp(4*x);

N = 5;

th = linspace(-pi, pi, N);

vals = f(th);

fft(vals)

% it has the same output as the julia fft!!