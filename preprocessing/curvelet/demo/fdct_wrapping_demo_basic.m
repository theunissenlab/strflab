disp(' ');
disp('fdct_wrapping_demo_basic.m -- Displays a curvelet')
disp ('both in the spatial and frequency domains.');
disp(' ');
disp(['This is achieved by setting all the coefficients in the curvelet'])
disp(['domain to zero except that at the required location (which'])
disp(['is set to one). The curvelet is obtained by taking the'])
disp(['adjoint curvelet transform. Notice how the curvelet is sharply '])
disp(['localized in both space and frequency.']); 
disp(' ');

% fdct_wrapping_demo_basic.m -- Displays a curvelet both in the spatial and frequency domains.

m = 32;
n = 32;

X = zeros(m,n);

%forward curvelet transform
disp('Take curvelet transform: fdct_wrapping');
tic; C = fdct_wrapping(X,0,1,3,8); toc;

%specify one curvelet
s = 3;
w = 1;
[A,B] = size(C{s}{w});
a = ceil((A+1)/2)+3;
b = ceil((B+1)/2)+4;
C{s}{w}(a,b) = 1;

%adjoint curvelet transform
disp('Take adjoint curvelet transform: ifdct_wrapping');
tic; Y = ifdct_wrapping(C,0,m,n); toc;

%display the curvelet
F = ifftshift(fft2(fftshift(Y)));
subplot(1,2,1); colormap gray; imagesc(real(Y)); axis('image'); ...
    title('a curvelet: spatial viewpoint');
subplot(1,2,2); colormap gray; imagesc(abs(F)); axis('image'); ...
    title('a curvelet: frequency viewpoint');

%get parameters
[SX,SY,FX,FY,NX,NY] = fdct_wrapping_param(C);

