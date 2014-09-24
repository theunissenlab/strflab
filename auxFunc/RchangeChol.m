function R = RchangeChol(R,type,j,strf,A,lambda,datIdx)

% function R = RchangeChol(R,type,j,strf,A,lambda,datIdx)
% Fast Cholesky insert and remove functions
%
% INPUT:
%		   R: current Cholesky matrix
%		type: 'insert' to update R in a Cholesky factorization 
%			  	R'R = X'X of a data matrix X
%			  'delete' to deletes a variable from the X'X matrix 
%				in a Cholesky factorisation R'R = X'X. Returns the 
%				downdated R. 
%		   j: index of variable that has just joined the active set
%		strf: strf structure, needed just for 'insert'
%		   A: active set, necessary just for 'insert'
%	  lambda: elastic net hyperparameter, necessary just for 'insert'
%	  datIdx: array of indices for the training set
%
% OUTPUT:
%			R: updated cholesky matrix
%
% SEE ALSO: trnLARS, computeLARSstep
%
% Adapted from the implementation of Karl Skoglund, IMM, DTU, kas@imm.dtu.dk
% to include delays 
% Lucas Pinto, December 2009, lucaspinto@berkeley.edu

global globDat

switch type
    case 'insert'
        if nargin < 7; datIdx = 1:size(globDat.stim,1); end
        x = globDat.stim(:,strf.varlookup(j,1)); 
        xdelay = strf.varlookup(j,2); 
		x = [zeros(xdelay,1);x(1:end-xdelay,1)]; 
        x = x(datIdx,:); % column vector representing the variable to be added  
        
        X = zeros(length(datIdx),length(A)); % stim matrix containing currently active variables
        for ti=1:length(A) % take care of the delays
            delay = strf.varlookup(A(ti),2);
            delaystim = [zeros(delay,1); globDat.stim(1:end-delay,strf.varlookup(A(ti),1))];
            X(:,ti) = delaystim(datIdx,:); %[zeros(delay,1); globDat.stim(1:datIdx(end)-delay,strf.varlookup(A(ti),1))];
        end
       
        diag_k = (x'*x + lambda)/(1 + lambda); % diagonal element k in X'X matrix
        if isempty(R)
            R = sqrt(diag_k);
        else
            col_k = x'*X/(1 + lambda); % elements of column k in X'X matrix
            R_k = R'\col_k'; % R'R_k = (X'X)_k, solve for R_k
            R_kk = sqrt(diag_k - R_k'*R_k); % norm(x'x) = norm(R'*R), find last element by exclusion
            R = [R R_k; [zeros(1,size(R,2)) R_kk]]; % update R
        end
        
    case 'delete'
        R(:,j) = []; % remove column j
        n = size(R,2);
        for k = j:n
            p = k:k+1;
            [G,R(p,k)] = planerot(R(p,k)); % remove extra element in column
            if k < n
                R(p,k+1:n) = G*R(p,k+1:n); % adjust rest of row
            end
        end
        R(end,:) = []; % remove zero'ed out row
end