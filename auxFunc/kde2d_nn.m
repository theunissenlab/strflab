function [pxy, bw] = kde2d_nn(data, xgrid, ygrid)

    K = 10;

    npts = size(data, 1);
    
    %% compute nearest-neighbor bandwidths for each point
    tic;
    %bw = compute_nn_bandwidths(data, K);
    bw = compute_nn_bandwidths2(data, K);
    nntime = toc;
    fprintf('Computing NN bandwidths took %f seconds\n', nntime);
    
    %% recenter bandwidth so it's >= 1
    %bw = (bw / max(bw(:))) + 1;
    %bw = exp(bw);
    %bw = bw / max(bw(:));
    %bw = log(bw);
    
    tic;
    %% evaluate pxy using NN Gaussian kernel
    [X, Y] = meshgrid(xgrid, ygrid);
    pxy = zeros(size(X));
    
    for k = 1:npts
        xcent = X - data(k, 1);
        ycent = Y - data(k, 2);
        bwx = bw(k, 1);
        bwy = bw(k, 2);
        
        h = 1;
        g = (xcent.^2 / (2*bwx)) + (ycent.^2 / (2*bwy));
        pxy = pxy + h*exp(-g);        
    end
    etime = toc;
    %fprintf('Evaluating NN took %f seconds\n', etime); 
    
    %% normalize p(x,y)
    pxy = pxy / sum(pxy(:));
    %{
    figure; hold on;
    imagesc(xgrid, ygrid, pxy); axis tight;
    plot(data(:, 1), data(:, 2), 'w.');
    colorbar;
    %}

end

function bw = compute_nn_bandwidths(data, K)

    bw = zeros(size(data));
    for k = 1:size(data, 1)
       
        pntx = data(k, 1);
        pnty = data(k, 2);      
        dx = (data(:, 1) - pntx).^2;
        dy = (data(:, 2) - pnty).^2;
        
        rdist = sqrt(dx+dy);
        
        % remove K-1 nearest neighbors (including [pntx, pnty])
        for m = 1:K
            [minVal, minIndx] = min(rdist);
            rdist(minIndx) = NaN;
        end
        
        % get top K nearest neighbor distances
        ndists = zeros(K, 2);
        for m = 1:K
            [minVal, minIndx] = min(rdist);            
            nbx = data(minIndx, 1);
            nby = data(minIndx, 2);
            ndists(m, :) = [nbx nby];
            rdist(minIndx) = NaN;
        end
        
        mbx = mean(ndists(:, 1));
        mby = mean(ndists(:, 2));
        
        % get vertical and horizontal components of distance vector
        bwx = abs(pntx - mbx);
        bwy = abs(pnty - mby);
        bw(k, :) = [bwx bwy];
    end    
end


function bw = compute_nn_bandwidths2(data, K)

    bw = zeros(size(data));    
    kdt = kdtree(data);
    
    x = data(:, 1);
    minx = min(x); maxx = max(x);
    xwidth = maxx - minx;    
    
    y = data(:, 2);
    miny = min(y); maxy = max(y);
    ywidth = maxy - miny;    
    
    xEps = 0.01*xwidth;
    yEps = 0.01*ywidth;
    
    npts = size(data, 1);
    
    dcnts = zeros(npts, 1);
        
    for k = 1:npts
    
        px = data(k, 1);
        py = data(k, 2);
        
        rng = [px - xEps, px + xEps; ...
               py - yEps, py + yEps];
        
        %% get neighbors within rectangle by searching kdtree
        cindxs = kdtree_range(kdt, rng);
        
        %% record # of neighboring points
        dcnts(k) = length(cindxs);
        
        %% compute distances
        dx = (data(cindxs, 1) - px).^2;
        dy = (data(cindxs, 2) - py).^2;        
        rdist = sqrt(dx+dy);
        
        % ignore (px, py)
        rdist(rdist == 0) = NaN;        
        % sort
        [rdist, rindx] = sort(rdist);
        %get top K       
        nk = min(K, length(rindx));
        if length(rindx) == 0
            error('[compute_nn_bandwidths2] No nearest neighbors found, increase the cube size!\n');
        end
        topk = [rv(dx(rindx(1:nk))) rv(dy(rindx(1:nk)))];
        
        %get means 
        mx = mean(topk(:, 1));
        my = mean(topk(:, 2));
        
        % get vertical and horizontal components of distance vector
        bwx = abs(px - mx);
        bwy = abs(py - my);
        bw(k, :) = [bwx bwy];        
    end    
    
    %% compute densities and sparseness
    dens = dcnts / (4*xEps*yEps);
    dens = dens / max(dens);
    
    sps = 1 ./ dens;
    sps = log(sps);
    sps = sps + max(abs(sps));
    
    %% weight bandwidths by sparseness
    %bw = bsxfun(@times, bw, sps);
        
end
