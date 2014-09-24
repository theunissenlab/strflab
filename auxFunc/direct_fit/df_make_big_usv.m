function [big_u,big_s,big_v,max_stimnorm,big_stim_mat,max_s] = df_make_big_usv(nb,nf,fstim)
%  This function returns a cell array of the U, S and V matrices generated
%  for a stimulus set.  These need to be calculated only once per dataset (and not per tol value),
%  so we can save time if we do them in a "do_cached_calc" call, even if
%  the dataset used is unique.  

clear j
ranktest=zeros(1,nf);
stimnorm=zeros(1,nf);
big_u = {};
big_s = {};
big_v = {};
big_stim_mat = {};
%disp('got here')
stim_mat = zeros(nb,nb)+j.*zeros(nb,nb);
% ========================================================
% Find the maximum norm of all the matrices
% ========================================================
for iff=1:nf
    nc = 1;
    for fb_indx=1:nb
        for fb_indx2=fb_indx:nb
            stim_mat(fb_indx,fb_indx2) = fstim(nc,iff);
            if fb_indx ~= fb_indx2
                stim_mat(fb_indx2,fb_indx) = conj(fstim(nc,iff));
            end
            nc = nc +1;
        end
    end
    stimnorm(1,iff)=norm(stim_mat);
end

% One way of doing it
%ranktol=tol*max(stimnorm);
max_stimnorm = max(stimnorm);
max_s = -999;
for iff=1:nf
    % Stuff stim matrix and cross-correlation vectors
    nc = 1;
    for fb_indx=1:nb
        for fb_indx2=fb_indx:nb
            stim_mat(fb_indx,fb_indx2) = fstim(nc,iff);
            if fb_indx ~= fb_indx2
                stim_mat(fb_indx2,fb_indx) = conj(fstim(nc,iff));
            end
            nc = nc +1;
        end
        % cross_vect(i) = conj(fstim_spike(i,iff));
%         cross_vect(fb_indx) = fstim_spike(fb_indx,iff);
%         for iJN=1:nJN
%             % =======================================================
%             %  JXZ: 7/12/2005
%             % cross_vectJN(iJN,i) = conj(stim_spike_JNf(i,iff,iJN));
%             % =======================================================
%             %cross_vectJN(iJN,fb_indx) = stim_spike_JNf(fb_indx,iff,iJN);
%             %Normalized by the duration of each individual stim
%             cross_vectJN(iJN,fb_indx) = (sum(durs)*fstim_spike(fb_indx,iff) - sum(durs([1:(iJN-1) (iJN+1):end]))*stim_spike_JNf(fb_indx,iff,iJN))/durs(iJN);
% 
%         end
    end

    % do an svd decomposition
    %ranktest(1,iff)=rank(stim_mat,ranktol);
    %[u,s,v] = svd(stim_mat);
    %[big_u{iff},big_s{iff},big_v{iff}] = do_cached_calc('svd',stim_mat);
    
    %  Use do_cash_calc here just for debugging; afterwards it will be more
    %  efficient to just calculate the svd since we don't expect individual
    %  svd results to be reusable.
    
    %[u,s,v] = do_cached_calc('svd',stim_mat);
    
    %  Not debugging anymore :-)
    [u,s,v] = svd(stim_mat);
    big_u{iff} = u;
    big_s{iff} = s;
    big_v{iff} = v;
    big_stim_mat{iff} = stim_mat;
    max_s = max(max_s,max(diag(s)));
    %[big_u{iff},big_s{iff},big_v{iff}] = svd(stim_mat);
end
%disp('even got here')
