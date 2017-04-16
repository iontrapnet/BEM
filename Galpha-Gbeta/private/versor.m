function [u,vmod] = versor(v)

% VERSOR calculates unit vectors of a set of vectors
% 
% USE:
% [u,vmod] = versor(v)
% 
% INPUT:
% 'v': set of vectors
% 
% OUTPUT:
% 'u': versors associated to v
% 'vmod': magnitute of vectors v
%
% NOTE:
%
% VERSION:
% Date: 03.03.2010
% Copyright(C) 2010-2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 19.03.2010: added exception handling when vmod == 0
% 16.09.2011: changes in help
% 14.09.2013: uses BSXFUN
% 18.09.2013: fixed bug when IDX is empty

u = zeros(size(v));
vmod = sqrt(sum(v.^2,2));
idx = find(vmod);
if ~isempty(idx)
    u(idx,:) = bsxfun(@rdivide,v(idx,:),vmod(idx));
end
