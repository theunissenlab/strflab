%%% 6 plots for 6 real basis functions
figure(1);clf;
subplot(321);
w1 = zeros(256); w2 = zeros(256);
w1(32,96) = 1+j;
x = icdwt2(w1,w2,2);
imagesc(abs(fftshift(fft2(x))));axis square;axis off;
line([94 162],[256 1]);
line([1 256],[94 162]);
title('+75');

subplot(322);
w1 = zeros(256); w2 = zeros(256);
w1(96,96) = 1;
x = icdwt2(w1,w2,2);
imagesc(abs(fftshift(fft2(x))));axis square;axis off;
line([1 256],[1 256]);
line([1 256],[256 1]);
title('+45');

subplot(323);
w1 = zeros(256); w2 = zeros(256);
w1(96,32) = 1;
x = icdwt2(w1,w2,2);
imagesc(abs(fftshift(fft2(x))));axis square;axis off;
line([94 162],[1 256]);
line([1 256],[162 94]);
title('+15');

subplot(324);
w1 = zeros(256); w2 = zeros(256);
w2(96,32) = 1;
x = icdwt2(w1,w2,2);
imagesc(abs(fftshift(fft2(x))));axis square;axis off;
line([94 162],[256 1]);
line([1 256],[94 162]);
title('-15');

subplot(325);
w1 = zeros(256); w2 = zeros(256);
w2(96,96) = 1;
x = icdwt2(w1,w2,2);
imagesc(abs(fftshift(fft2(x))));axis square;axis off;
line([1 256],[1 256]);
line([1 256],[256 1]);
title('-45');

subplot(326);
w1 = zeros(256); w2 = zeros(256);
w2(32,96) = 1;
x = icdwt2(w1,w2,2);
imagesc(abs(fftshift(fft2(x))));axis square;axis off;
line([94 162],[1 256]);
line([1 256],[162 94]);
title('-75');

%%% 1 basis function with 6 orientation lines
figure(1);clf;
w1 = zeros(256); w2 = zeros(256);
w1(32,96) = 1;
x = icdwt2(w1,w2,2);
imagesc(abs(fftshift(fft2(x))));axis square;axis off;
line([94 162],[256 1]);
line([1 256],[94 162]);
line([1 256],[1 256]);
line([1 256],[256 1]);
line([94 162],[1 256]);
line([1 256],[162 94]);
title('+75');

%%% Linear combinations of basis functions
w1 = zeros(256); w2 = zeros(256);
w1(32,96) = 1;
w1(96,96) = 1;
w1(96,32) = 1;
w2(32,96) = 1;
w2(96,96) = 1;
w2(96,32) = 1;
x = icdwt2(w1,w2,2);
imagesc(abs(fftshift(fft2(x))));axis square;axis off;
imagesc(angle(fftshift(fft2(x))));axis square;axis off;
temp = angle(fftshift(fft2(x)));
imagesc(mod(temp,pi));axis square;axis off;

%%% Basis functions at different scales
w1 = zeros(256); w2 = zeros(256);
w1(32,96) = 1;w1(96,96) = 1;w1(96,32) = 1;
w2(32,96) = 1;w2(96,96) = 1;w2(96,32) = 1;
w3 = zeros(256); w4 = zeros(256);
w3(16,48) = 1;w3(48,48) = 1;w3(48,16) = 1;
w4(16,48) = 1;w4(48,48) = 1;w4(48,16) = 1;
x1 = icdwt2(w1,w2,2);
x3 = icdwt2(w3,w4,3);
figure(1);clf;
imagesc(abs(fftshift(fft2(x1))));axis square;axis off;
line([65 192],[65 65]);
line([65 65],[65 192]);
line([65 192],[192 192]);
line([192 192],[65 192]);
line([97 160],[97 97]);
line([97 97],[97 160]);
line([97 160],[160 160]);
line([160 160],[97 160]);
figure(2);clf;
imagesc(abs(fftshift(fft2(x3))));axis square;axis off;
line([97 160],[97 97]);
line([97 97],[97 160]);
line([97 160],[160 160]);
line([160 160],[97 160]);
line([113 144],[113 113]);
line([113 113],[113 144]);
line([113 144],[144 144]);
line([144 144],[113 144]);

%%% Testing ratios of frequency response, for fixed angle
w1 = zeros(256); w2 = zeros(256);
w1(32,96) = 1 + j;
xp75 = fftshift(fft2(icdwt2(w1,w2,2)));
w1 = zeros(256); w2 = zeros(256);
w1(96,96) = 1 + j;
xp45 = fftshift(fft2(icdwt2(w1,w2,2)));
w1 = zeros(256); w2 = zeros(256);
w1(96,32) = 1 + j;
xp15 = fftshift(fft2(icdwt2(w1,w2,2)));
w1 = zeros(256); w2 = zeros(256);
w2(96,32) = 1 + j;
xm15 = fftshift(fft2(icdwt2(w1,w2,2)));
w1 = zeros(256); w2 = zeros(256);
w2(96,96) = 1 + j;
xm45 = fftshift(fft2(icdwt2(w1,w2,2)));
w1 = zeros(256); w2 = zeros(256);
w2(32,96) = 1 + j;
xm75 = fftshift(fft2(icdwt2(w1,w2,2)));
figure(1);clf;
l75 = zeros(1,128);
l45 = zeros(1,128);
l15 = zeros(1,128);
for ii = 1:128
  l75(ii) = xp75(ii,ii);
  l45(ii) = xp45(ii,ii);
  l15(ii) = xp15(ii,ii);
end
plot([1:128],abs(l75),[1:128],abs(l45));
plot([1:128],angle(l75),[1:128],angle(l45));
legend('+75','+45');
title('Cross section of Fourier transform, angle +45 degrees');
figure(2);clf;
l75 = zeros(1,64);
l45 = zeros(1,64);
l15 = zeros(1,64);
for ii = 1:64
  l75(ii) = xp75(2*ii,ii+64);
  l45(ii) = xp45(2*ii,ii+64);
  l15(ii) = xp15(2*ii,ii+64);
end
plot([1:64],abs(l75),[1:64],abs(l45),[1:64],abs(l15));
plot([1:64],angle(l75),[1:64],angle(l45),[1:64],angle(l15));
legend('+75','+45','+15');
title('Cross section of Fourier transform, angle atan(1/2) degrees');



