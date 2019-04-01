%% let's see if I know what's going on with Fourier stuff...

N = 101;
a = 1/2;
b = 1/8;
x0 = 1;
x1 = -1/16;

x = linspace(-1, 1, N);
f = zeros(1,N);
f(abs(x) < a) = x0;
f(abs(x) < b) = f(abs(x) < b) + x1;

figure; plot(x, f);

%% Fourier coefficients
k = 100;
fk = x0./(pi*(1:k)) .* sin(pi * a * (1:k)) + ...
    x1./(pi*(1:k)) .* sin(pi * b * (1:k));

figure; plot(1:k, fk);

% construct a k-by-N matrix representing basis functions
K = repmat((1:k)', 1, N);
basis = exp(1i * pi * K .* x);
fhat = fk * basis;

figure; plot(x, real(fhat));