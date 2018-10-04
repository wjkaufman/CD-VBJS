function [row_a, col_a, row_d, col_d] = glrt_cd(image1, image2, tz, w_size)

% Developed by Jason Mossing AFRL/RYAS at WPAFB 
% Distribution Statement A: Approved for public release; distribution is 
% unlimited. 88 ABW-10-4197 
%Equation based of off Leslie Novak's equation 25 in "Change Detection for 
% multi-polarization, multi-pass SAR" proceedings of SPIE (2005)
% http://lesnovak.com/images/SPIE05_Novak.pdf 
%Inputs:
% image1 - This is a complex image (2-D matrix) is the reference image
% image2 - This is a complex image (2-D matrix) is the test image
% tz - This is the threshold
% w_size - This is length or width of the window...w_size squared is the
% number of pixels in window
%Outputs:
% row_a - This is the row of the arrivals
% col_a - This is the column of the arrivals
% row_d - This is the row of the departures
% col_d - This is the column of the departures
%Description:
% Performs Generalized Likelihood Ratio Test, Plots the arrivals (blue) and
% departures (red) on top of the referance image

w=ones(w_size,w_size);
% creating G which contains values to be compared to the threshold
c1=conv2(abs(image1).^2, w, 'same');

c2=conv2(abs(image2).^2, w, 'same');

clear image2

G=c2./c1;

clear c1 c2 w w_size

% finding arrivals

[row_a, col_a] = find(G > tz); 
row_a=single(row_a); 
col_a=single(col_a);

% finding departures
[row_d, col_d] = find(G < (1/tz)); 
row_d=single(row_d); 
col_d=single(col_d);

clear G tz

% show results

imagesc(20*log10(abs(image1)));
% used for visual aid
colormap(bone);
axis square xy;
colorbar
ax = caxis - 10;
ax(1) = ax(2) - 70;
caxis(ax);
% plot arrivals and departures, red fled blue new
hold on
plot(col_d,row_d,'r.','markersize',0.01);
plot(col_a,row_a,'b.','markersize',0.01);
hold off