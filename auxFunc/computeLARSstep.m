function [gamma, u1, w] = computeLARSstep(grad,R,A,I,C,d1,d2,datIdx,strf)

% function [gamma, u1, w] = computeLARSstep(grad,R,A,I,C,d1,d2,datIdx,strf)
% computes step for the LARS algorithm. Called by trnLARS
%
% INPUT:
% 		grad: gradient
%		   R: Choletsky matrix
%		   A: Active set 
%		   I: Inactive set	
%		   C: maximal correlation
%		  d1: hyperparameter for elastic net
%		  d2: hyperparameter for elastic net
%	  datIdx: array of indices for the training set
%		strf: strf structure
%
% OUTPUT:
%		[gamma]: LARS step
%		[u1, w]: variables used to calculate equiangular vector  
%
% SEE ALSO: trnLARS
%
% Adapted from the implementation of Karl Skoglund, IMM, DTU, kas@imm.dtu.dk
% to include delays 
% Lucas Pinto, December 2009, lucaspinto@berkeley.edu

global globDat

vars = length(A);
nvars = strf.nWts -1;
ndelays = length(strf.delays);
t = length(datIdx);

s = sign(grad(A))'; % get the signs of the correlations
GA1 = R\(R'\s);
AA = 1/sqrt(sum(GA1.*s));
w = AA*GA1;

stim = zeros(t,vars);
for ti=1:vars 
    delay = strf.varlookup(A(ti),2); %keyboard
    delaystim = [zeros(delay,1); globDat.stim(1:end-delay,strf.varlookup(A(ti),1))];
    stim(:,ti) = delaystim(datIdx,:);%[zeros(delay,1); globDat.stim(1:datIdx(end-delay),strf.varlookup(A(ti),1))];
end

u1 = [stim*w*d2; zeros(strf.delays(end),1)]; % equiangular direction (unit vector) part 1
% u1 = stim*w*d2;

u2 = zeros(nvars, 1); u2(A) = d1*d2*w; % part 2
if vars == nvars % if all variables active, go all the way to the lsq solution
    gamma = C/AA;
else
    u1m=zeros(nvars,1);
    for iii = 1:ndelays
        % tic
        % delaystim = [zeros(strf.delays(iii),size(globDat.stim,2)); globDat.stim(1:end-strf.delays(iii),:)];
        % a=toc
        % u1m(strf.nIn*(iii-1)+1:strf.nIn*iii) = delaystim(datIdx,:)'*u1;
        % keyboard
        % tic
        u1m(strf.nIn*(iii-1)+1:strf.nIn*iii) = globDat.stim(datIdx,:)'*u1(1+strf.delays(iii):t+strf.delays(iii));
        % b=tocmax
        % keyboard
            %[zeros(strf.delays(iii),size(globDat.stim,2)); globDat.stim(1:datIdx(end)-strf.delays(iii),:)]'*u1; 
    end

    a = (u1m + d1*u2)*d2; % correlation between each variable and equiangular vector
    temp = [(C - grad(I)')./(AA - a(I)); (C + grad(I)')./(AA + a(I))];
    gamma = min([temp(temp > 0); C/AA]);
end
u1 = u1(1:t);