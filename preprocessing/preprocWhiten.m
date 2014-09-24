function [wstim, params, winfo] = preprocWhiten(stim, params)

    if ~isfield(params, 'method')
        params.method = 'meansub';
    end

    cmat = cov(stim);
    
    [V, D] = eig(cmat);
    
    Dsqrt = diag(diag(D).^ -0.5);
    
    switch params.method
        
        case 'meansub'
            stim = bsxfun(@minus, stim, mean(stim, 1));
            tmat = Dsqrt * V';            
            
        case 'matsq'
            tmat = V*Dsqrt*V';       
    end
    
    wstim = zeros(size(stim));
    for k = 1:size(stim, 1)
        wstim(k, :) = (tmat*stim(k, :)')';
    end

    wmean = mean(wstim, 1);
    wstim = bsxfun(@minus, wstim, wmean);
    
    %{
    wcmat = cov(wstim);
    
    figure; hold on;
    subplot(2, 1, 1); hold on;
    imagesc(cmat); axis tight;
    title('Covariance Matrix');
    colorbar;
    
    subplot(2, 1, 2); hold on;
    imagesc(wcmat); axis tight;
    title('Whitened Covariance Matrix');
    colorbar;
    %}
    
    winfo = struct;
    winfo.mean = wmean;
    winfo.cmat = cmat;
    