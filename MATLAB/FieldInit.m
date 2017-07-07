function [cb,xr,yr,zr,pb] = FieldInit(path, varargin)
    %global cb xr yr zr pb
    cb = ChargeBasis(path);
    if nargin == 4
        xr = varargin{1};
        yr = varargin{2};
        zr = varargin{3};
        pb = PotentialBasis(path, xr, yr, zr);
    else
        pb = 0;
    end
    %ans = pb;
end