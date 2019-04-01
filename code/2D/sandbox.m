N=100;
x = linspace(-2*pi,2*pi,N);
y = linspace(-2*pi,2*pi,N)';

f= sin(3*x+.4) + cos(5*y) + .2;
figure; imagesc(f);

Y = fft2(f)/sqrt(numel(f));

figure; imagesc(real(Y));
figure; plot(real(Y(:,1)));