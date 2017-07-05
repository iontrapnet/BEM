function varargout = EqualRangeInterpolation(point, xr, yr, zr, varargin)
    n = [xr(3),yr(3),zr(3)]+1;
    d = [];
    if n(1) > 1
        d(1) = (n(1)-1)*(point(1)-xr(1))/(xr(2)-xr(1));
    else
        d(1) = 0;
    end
    if n(2) > 1
        d(2) = (n(2)-1)*(point(2)-yr(1))/(yr(2)-yr(1));
    else
        d(2) = 0;
    end
    if n(3) > 1
        d(3) = (n(3)-1)*(point(3)-zr(1))/(zr(2)-zr(1));
    else
        d(3) = 0;
    end
    i = floor(d);
    d = d - i;
    coe = [(1-d(1))*(1-d(2))*(1-d(3)) (1-d(1))*(1-d(2))*d(3) (1-d(1))*d(2)*(1-d(3)) (1-d(1))*d(2)*d(3) d(1)*(1-d(2))*(1-d(3)) d(1)*(1-d(2))*d(3) d(1)*d(2)*(1-d(3)) d(1)*d(2)*d(3)];
    idx = 1+[i(3)+i(2)*n(3)+i(1)*n(2)*n(3) 1+i(3)+i(2)*n(3)+i(1)*n(2)*n(3) i(3)+(1+i(2))*n(3)+i(1)*n(2)*n(3) 1+i(3)+(1+i(2))*n(3)+i(1)*n(2)*n(3) i(3)+i(2)*n(3)+(1+i(1))*n(2)*n(3) 1+i(3)+i(2)*n(3)+(1+i(1))*n(2)*n(3) i(3)+(1+i(2))*n(3)+(1+i(1))*n(2)*n(3) 1+i(3)+(1+i(2))*n(3)+(1+i(1))*n(2)*n(3)];
    for k=1:(nargin-4)
        varargout{k} = coe*varargin{k}(idx)';
    end
end