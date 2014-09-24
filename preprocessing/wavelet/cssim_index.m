function [cssim, band_cssim] = cssim_index(im1, im2, level, or, guardb)

[pyr1, pind] = buildSCFpyr(im1, level, or-1);
[pyr2, pind] = buildSCFpyr(im2, level, or-1);


%window = fspecial('gaussian', [11 11], 2);
window = ones(7);
window = window./sum(sum(window));
C = 0;
gb = guardb/(2^(level-1));

for i=1:or
   bandind = i+(level-1)*or+1;
   band1 = pyrBand(pyr1, pind, bandind);
   band2 = pyrBand(pyr2, pind, bandind);
   band1 = band1(gb+1:end-gb, gb+1:end-gb);
   band2 = band2(gb+1:end-gb, gb+1:end-gb);
   
   corr = band1.*conj(band2);
   varr = abs(band1).^2 + abs(band2).^2;
   corr_band = filter2(window, corr, 'valid');
   varr_band = filter2(window, varr, 'valid');
   cssim_map = (2*abs(corr_band) + C)./(varr_band + C);

   band_cssim(i) = mean2(cssim_map);
end

cssim = mean(band_cssim);
