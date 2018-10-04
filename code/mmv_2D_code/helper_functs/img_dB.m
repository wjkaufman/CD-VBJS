function img_out = imgDb(img_in)
% input img_in in the size you want img_out to be (NxM) 

img_out = 20*log10(abs(img_in)/max(max(abs(img_in)))); 

end