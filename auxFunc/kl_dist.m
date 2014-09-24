function d = kl_dist(px, py, rng)

    thresh = 0;
    nz = px > thresh;
    
    dx = rng(2)-rng(1);
    
    peps = 1e-100;
    pxnz = px(nz);
    pynz = py(nz) + peps;
        
    %lvals = pxnz.*(log2(pxnz) - log2(pynz));
    %lvals = pxnz.*(log2(pxnz ./ pynz));
    %d1 = sum(lvals)*dx
    
    lvals1 = pxnz.*log2(pxnz);
    lvals2 = pxnz.*log2(pynz);
    d = (sum(lvals1) - sum(lvals2))*dx;
    