clear all;
close all;

% Read the original image
load lenna
im=X;
[M N] = size(im);

% Create the first image to be compared
im1 = 128*ones(M+32, N+32);
im1(16+1:16+M, 16+1:16+N) = im;

% Create the second image to be compared
% Distortions include spatical shift, contrast change, luminance change,
% and Gaussian noise contamination
sx = -3 + 3;
sy = -2 - 1;
Mr = M+6;
Nr = N+4;
imr = imresize(im, [Mr Nr]);
im2 = 128*ones(M+32, N+32);
im2(16+1+sy:16+Mr+sy, 16+1+sx:16+Nr+sx) = 1.3*(imr-128)+128+30;
im2 = im2 + 10*randn(size(im2));
%g = fspecial('gaussian', [11 11], 1.2);
%im2 = filter2(g, im2, 'same');

figure;
imshow(im1/255);
figure;
imshow(im2/255);

% Compare the images using MSE and standard SSIM 
gb = 32;
img1 = im1(gb+1:gb+M, gb+1:gb+N);
img2 = im2(gb+1:gb+M, gb+1:gb+N);
mse = mean2((img1 - img2).^2)
ssim = ssim_index(img1, img2)

% CW-SSIM comparison at level 1 (finest scale)
level = 1;
or = 4;
[cssim band_cssim] = cssim_index(im1, im2, level, or, gb);
cssim
band_cssim

% CW-SSIM comparison at level 2 (second finest scale)
level = 2;
or = 4;
[cssim band_cssim] = cssim_index(im1, im2, level, or, gb);
cssim
band_cssim

% CW-SSIM comparison at level 3 (third finest scale)
level = 3;
or = 4;
[cssim band_cssim] = cssim_index(im1, im2, level, or, gb);
cssim
band_cssim