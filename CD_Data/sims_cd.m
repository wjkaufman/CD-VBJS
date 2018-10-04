function [row_a, col_a, row_d, col_d] = sims_cd(image1, image2, tz)

% Simple Image Magnitutde Subtraction
% Distribution Statement A: Approved for public release; distribution is 
% unlimited. 88 ABW-10-4197
% Developed by Jason Mossing AFRL/RYAS at WPAFB
%Inputs:
% image1 - This is a complex image (2-D matrix) is the refrence image
% image2 - This is a complex image (2-D matrix) is the test image
%Outputs:
% row_a - This is the row of the arrivals
% col_a - This is the column of the arrivals
% row_d - This is the row of the departures
% col_d - This is the column of the departures
%Description:
% Subtracts the absolute values of image1 from image2 and if over tz it is
% an arrival and if under -tz it is a departure, and plots arrivals and
% departures over the refrence image

% find a which contains values to be compared to the threshold
a = abs(image2) - abs(image1);
% imagesc refrence image
imagesc(20*log10(abs(image1)));
% used for visual aid
colormap(bone);
axis square xy;
colorbar
ax = caxis - 10;
ax(1) = ax(2) - 70;
caxis(ax);
hold on
% find and plot arivals
[row_a,col_a]=find(a >= tz);
plot(col_a,row_a,'b.','markersize',1);
% find and plot departures
[row_d,col_d]=find(a <= -tz);
plot(col_d,row_d,'r.','markersize',1);