function toeMatrix = calcSparseToeplitz(stim,sampleIdxs,extentForwards,extentBackwards)
    sampleTimesSel = sampleIdxs( (sampleIdxs > extentForwards) & (sampleIdxs < length(stim)-extentBackwards+1));
    toeMatrix = zeros(length(sampleTimesSel),extentBackwards+extentForwards+1);
    for ii=1:length(sampleTimesSel)
        % note, this is just an explicit computation of the cross-correlation
        toeMatrix(ii,:) = stim(sampleTimesSel(ii)+[-extentForwards:extentBackwards]);
    end
    toeMatrix = fliplr(toeMatrix);
end