N=64;
x = linspace(-2*pi,2*pi,N);
y = linspace(-2*pi,2*pi,N)';

f= sin(5*x) + cos(20*y) + 5;
figure; imagesc(f);

Y = fft2(f)/sqrt(numel(f));

figure; imagesc(real(Y));