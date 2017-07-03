function varargout = EqualRangeInterpolation(point, xr, yr, zr, varargin)
    nx = xr(3) + 1;
    ny = yr(3) + 1;
    nz = zr(3) + 1;
    d = [];
    if nx > 1
        d(1) = (nx-1)*(point(1)-xr(1))/(xr(2)-xr(1));
    else
        d(1) = 0;
    end
    if ny > 1
        d(2) = (ny-1)*(point(2)-yr(1))/(yr(2)-yr(1));
    else
        d(2) = 0;
    end
    if nz > 1
        d(3) = (nz-1)*(point(3)-zr(1))/(zr(2)-zr(1));
    else
        d(3) = 0;
    end
    i = floor(d);
    d = d - i;
    coe = [(1-d(1))*(1-d(2))*(1-d(3)) (1-d(1))*(1-d(2))*d(3) (1-d(1))*d(2)*(1-d(3)) (1-d(1))*d(2)*d(3) d(1)*(1-d(2))*(1-d(3)) d(1)*(1-d(2))*d(3) d(1)*d(2)*(1-d(3)) d(1)*d(2)*d(3)];
    idx = 1+[i(3)+i(2)*nz+i(1)*ny*nz 1+i(3)+i(2)*nz+i(1)*ny*nz i(3)+(1+i(2))*nz+i(1)*ny*nz 1+i(3)+(1+i(2))*nz+i(1)*ny*nz i(3)+i(2)*nz+(1+i(1))*ny*nz 1+i(3)+i(2)*nz+(1+i(1))*ny*nz i(3)+(1+i(2))*nz+(1+i(1))*ny*nz 1+i(3)+(1+i(2))*nz+(1+i(1))*ny*nz];
    for k=1:(nargin-4)
        varargout{k} = coe*varargin{k}(idx)';
    end
end