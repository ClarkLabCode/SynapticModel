function [ stimArray ] = MovingEdges(params, barParam, v)
%A single moving edge

pdEdge = 2*double((params.mask>0) & (v*params.t - params.x) > 0)-1;

ndEdge = fliplr(pdEdge);

stimArray = barParam.mlum + barParam.c * double(cat(3, pdEdge, -pdEdge, ndEdge, -ndEdge));

end

