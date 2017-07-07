function [pb,cb] = FieldInit(path, varargin)
    if nargin < 4
        [cb,pb] = ChargeBasis(path);
    else
        xr = varargin{1};
        yr = varargin{2};
        zr = varargin{3};
        pb = PotentialBasis(path, xr, yr, zr);
        cb = [xr;yr;zr];
    end
end