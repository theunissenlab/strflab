Dual-Tree Complex Wavelet Transform Pack - version 2.1

Nick Kingsbury and Cian Shaffrey, Cambridge University, September 2000.

Version 2.0 files (dated August 2000) have been optimised for speed
and memory use, principally by performing column filtering with CONV2 
instead of row filtering (which is significantly slower) and by avoiding 
unneccessary copying of matrices within functions.  The ROWFILT functions
have been replaced by COLFILT equivalents, and some lower level functions
have been removed.

Version 2.1 files (dated Sept 2000) contain 2 new functions, CWTBAND6 and
ICWTBAND6, and some bugs in the comments of other files have been corrected.
These new functions provide a more convenient way to handle complete sets
of 6 complex bandpass subbands from a given level of the transform.

The main functions:

dtwavedec => 1D DTCWT decomposition
dtwaverec => 1D DTCWT reconstruction
cwtband   => Extracts individual subbands from the result of DTWAVEDEC 
             at specified levels
icwtband  => Allows insertion of individual subbands into the input vector
             for DTWAVEREC.

dtwavedec2 => 2D DTCWT decomposition
dtwaverec2 => 2D DTCWT reconstruction
cwtband2   => Extracts individual subimages from the result of DTWAVEDEC2
              at specified levels
icwtband2  => Allows insertion of individual subimages into the input vector
              for DTWAVEREC2
cwtband6   => Extracts a set of 6 subimages from the result of DTWAVEDEC2
              at specified levels
icwtband6  => Allows insertion of a set of 6 subimages into the input vector
              for DTWAVEREC2.

Lower level functions:

colfilt   => Column filtering of a matrix with symmetric extension 
coldfilt  => Column filtering with decimation by 2.
colifilt  => Column filtering with interpolation (upsampling) by 2.
reflect   => Reflect a vector about max and min limits (used for sym extension).
draw      => Draw an image in a correctly sized figure window.
cimage5   => Draw a complex subimage using a colour palette for the complex numbers.

Various .MAT files contain the complex wavelet filter coefficients.


% To test the 1D case:

X = rand(512,1);
[C,L] = dtwavedec(X,3,'near_sym_b','qshift_b');
Z = dtwaverec(C,L,'near_sym_b','qshift_b');
max(abs(Z(:)-X(:))) % Test for approx zero error.

% And for vectors that are NOT divisable by 4:

X = [zeros(225,1); ones(225,1)];
[C,L] = dtwavedec(X,3,'near_sym_b','qshift_b');
Z = dtwaverec(C,L,'near_sym_b','qshift_b');
max(abs(Z(:)-X(:)))

l = cwtband(C,L,3,'l');
h1 = cwtband(C,L,1,'h');
h2 = cwtband(C,L,2,'h','real');

% To insert some modified subband (eg h1*2) back into C

[cm,V] = icwtband(h1*2,L,1,'h');
C(V) = cm;



% To test the 2D case:

load lenna   % This gives X as 256x256.
[C,S] = dtwavedec2(X,4,'antonini','qshift_c');
Z = dtwaverec2(C,S,'antonini','qshift_c');
max(abs(Z(:)-X(:))) % Test for approx zero error.

% Testing non-uniform sections:

[C,S] = dtwavedec2(X(1:120,1:130),4,'antonini','qshift_c');
Z = dtwaverec2(C,S,'antonini','qshift_c');
max(max(abs(Z-X(1:120,1:130))))

% Pick out separate complex subbands, display them using a complex
% colour palette, and reconstruct the image gradually,
% starting with the coarsest levels:

load lenna
[C,S] = dtwavedec2(X,4,'antonini','qshift_c');

% Pick out 6 subbands at each level and display them.
b1 = cwtband6(C,S,1);
figure; cimage5([b1(:,:,1) b1(:,:,3) b1(:,:,5);b1(:,:,2) b1(:,:,4) b1(:,:,6)])

b2 = cwtband6(C,S,2);
figure; cimage5([b2(:,:,1) b2(:,:,3) b2(:,:,5);b2(:,:,2) b2(:,:,4) b2(:,:,6)])

b3 = cwtband6(C,S,3);
figure; cimage5([b3(:,:,1) b3(:,:,3) b3(:,:,5);b3(:,:,2) b3(:,:,4) b3(:,:,6)])

b4 = cwtband6(C,S,4);
figure; cimage5([b4(:,:,1) b4(:,:,3) b4(:,:,5);b4(:,:,2) b4(:,:,4) b4(:,:,6)])

l4 = cwtband2(C,S,4,'l','real'); 
figure; draw(l4); drawnow       % Display final 'real' lowpass image in monochrome.

% Now reconstruct in D and display after each new level of subbands is added in.
D = C * 0;
[cm,V] = icwtband2(l4,S,4,'l','real'); D(V)=cm; % Level 4 lo-lo band.
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow

[cm,V] = icwtband6(b4,S,4);D(V) = cm;
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow

[cm,V] = icwtband6(b3,S,3);D(V) = cm;
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow

[cm,V] = icwtband6(b2,S,2);D(V) = cm;
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow

[cm,V] = icwtband6(b1,S,1);D(V) = cm;
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow
max(abs(Z(:)-X(:))) % Test for approx zero error.

% ***********************************************************************
% The following code is an older and more complicated way to select and
% insert subbands using just the cwtband2 and icwtband2 functions.

% Pick out separate complex subbands, display them using a complex
% colour palette, and reconstruct the image gradually,
% starting with the coarsest levels:

load lenna
[C,S] = dtwavedec2(X,4,'antonini','qshift_c');

% Pick out subbands.
h1 = cwtband2(C,S,1,'h','cplx');
v1 = cwtband2(C,S,1,'v','cplx');
d1 = cwtband2(C,S,1,'d','cplx');
figure; cimage5([v1 d1 h1]); drawnow   % Display level 1 complex subbands in colour.

h2 = cwtband2(C,S,2,'h','cplx');
v2 = cwtband2(C,S,2,'v','cplx');
d2 = cwtband2(C,S,2,'d','cplx');
figure; cimage5([v2 d2 h2]); drawnow   % Level 2

h3 = cwtband2(C,S,3,'h','cplx');
v3 = cwtband2(C,S,3,'v','cplx');
d3 = cwtband2(C,S,3,'d','cplx');
figure; cimage5([v3 d3 h3]); drawnow   % Level 3

h4 = cwtband2(C,S,4,'h','cplx');
v4 = cwtband2(C,S,4,'v','cplx');
d4 = cwtband2(C,S,4,'d','cplx');
figure; cimage5([v4 d4 h4]); drawnow   % Level 4

l4 = cwtband2(C,S,4,'l','real'); 
figure; draw(l4); drawnow              % Display final 'real' lowpass image in monochrome.

% Now reconstruct in D and display after each new level of subbands is added in.
D = C * 0;
[cm,V] = icwtband2(l4,S,4,'l','real'); D(V)=cm; % Level 4 lo-lo band.
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow

[cm,V] = icwtband2(h4,S,4,'h','cplx'); D(V)=cm; % Level 4 hi bands.
[cm,V] = icwtband2(v4,S,4,'v','cplx'); D(V)=cm;
[cm,V] = icwtband2(d4,S,4,'d','cplx'); D(V)=cm;
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow

[cm,V] = icwtband2(h3,S,3,'h','cplx'); D(V)=cm; % Level 3 hi bands.
[cm,V] = icwtband2(v3,S,3,'v','cplx'); D(V)=cm;
[cm,V] = icwtband2(d3,S,3,'d','cplx'); D(V)=cm;
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow

[cm,V] = icwtband2(h2,S,2,'h','cplx'); D(V)=cm; % Level 2 hi bands.
[cm,V] = icwtband2(v2,S,2,'v','cplx'); D(V)=cm;
[cm,V] = icwtband2(d2,S,2,'d','cplx'); D(V)=cm;
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow

[cm,V] = icwtband2(h1,S,1,'h','cplx'); D(V)=cm; % Level 1 hi bands.
[cm,V] = icwtband2(v1,S,1,'v','cplx'); D(V)=cm;
[cm,V] = icwtband2(d1,S,1,'d','cplx'); D(V)=cm;
Z = dtwaverec2(D,S,'antonini','qshift_c'); figure; draw(Z); drawnow
max(abs(Z(:)-X(:))) % Test for approx zero error.

********************************

Further information on the DT CWT can be obtained from papers
downloadable from NGK's website (given below). The best tutorial is in
the 1999 Royal Society Paper. In particular this explains the conversion
between 'real' quad-number subimages and pairs of complex subimages, which
is carried out in CWTBAND2 and ICWTBAND2 if 'cplx' is selected as the
required subimage mode. The Q-shift filters are explained in the ICIP
2000 paper and the paper for the Journal on Applied Computation and
Harmonic Analysis.

Cian Shaffrey and Nick Kingsbury, 
Cambridge University, August 2000

***********************************************************
Dr N G Kingsbury,
  Dept. of Engineering, University of Cambridge,
    Trumpington St., Cambridge CB2 1PZ, UK.
                      or
    Trinity College, Cambridge CB2 1TQ, UK.
Phone: (0 or +44) 1223 338514 / 332647;  Home: 1954 211152;
Fax: 1223 338564 / 332662;  E-mail: ngk@eng.cam.ac.uk
Web home page: http://www.eng.cam.ac.uk/~ngk/
***********************************************************

