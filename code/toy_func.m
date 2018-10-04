function fk = toy_func(x0, a, k)
% returns Fourier coefficients corresponding to toy func
% x0: height of hat
% a: width of hat
% k: wavenumber

fk = 2*x0./(pi*(k)) .* sin(pi * a * (k));
fk(k==0) = a*x0;

end