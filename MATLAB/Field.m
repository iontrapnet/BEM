function fields = Field(voltages, points, pb, cb)
    voltages = voltages';
    if size(pb, 1) < 4
        [~, ex, ey, ez] = Potential(pb, cb*voltages, points);
        fields = [ex;ey;ez]';
    else
        nx = cb(1,3)+1;
        ny = cb(2,3)+1;
        nz = cb(3,3)+1;
        noe = size(pb, 3);
        nop = size(points, 1);
        %d = [(nx-1)*(points(:,1)-xr(1))/(xr(2)-xr(1)) (ny-1)*(points(:,2)-yr(1))/(yr(2)-yr(1)) (nz-1)*(points(:,3)-zr(1))/(zr(2)-zr(1))]; 
        %i = floor(d);
        %d = d - i;
        %coe = [(1-d(:,1)).*(1-d(:,2)).*(1-d(:,3)) (1-d(:,1)).*(1-d(:,2)).*d(:,3) (1-d(:,1)).*d(:,2).*(1-d(:,3)) (1-d(:,1)).*d(:,2).*d(:,3) d(:,1).*(1-d(:,2)).*(1-d(:,3)) d(:,1).*(1-d(:,2)).*d(:,3) d(:,1).*d(:,2).*(1-d(:,3)) d(:,1).*d(:,2).*d(:,3)];
        %idx = 1+[i(:,3)+i(:,2)*nz+i(:,1)*ny*nz 1+i(:,3)+i(:,2)*nz+i(:,1)*ny*nz i(:,3)+(1+i(:,2))*nz+i(:,1)*ny*nz 1+i(:,3)+(1+i(:,2))*nz+i(:,1)*ny*nz i(:,3)+i(:,2)*nz+(1+i(:,1))*ny*nz 1+i(:,3)+i(:,2)*nz+(1+i(:,1))*ny*nz i(:,3)+(1+i(:,2))*nz+(1+i(:,1))*ny*nz 1+i(:,3)+(1+i(:,2))*nz+(1+i(:,1))*ny*nz];
        %fields = reshape(reshape(pb(2:4, idx, :),[24*nop noe])*voltages,[3*nop 8])*coe';
        fields = zeros(nop, 3);
        for k=1:nop
            d = [(nx-1)*(points(k,1)-cb(1,1))/(cb(1,2)-cb(1,1)) (ny-1)*(points(k,2)-cb(2,1))/(cb(2,2)-cb(2,1)) (nz-1)*(points(k,3)-cb(3,1))/(cb(3,2)-cb(3,1))]; 
            i = floor(d);
            d = d - i;
            coe = [(1-d(1))*(1-d(2))*(1-d(3));(1-d(1))*(1-d(2))*d(3);(1-d(1))*d(2)*(1-d(3));(1-d(1))*d(2)*d(3);d(1)*(1-d(2))*(1-d(3));d(1)*(1-d(2))*d(3);d(1)*d(2)*(1-d(3));d(1)*d(2)*d(3)];
            idx = 1+[i(3)+i(2)*nz+i(1)*ny*nz 1+i(3)+i(2)*nz+i(1)*ny*nz i(3)+(1+i(2))*nz+i(1)*ny*nz 1+i(3)+(1+i(2))*nz+i(1)*ny*nz i(3)+i(2)*nz+(1+i(1))*ny*nz 1+i(3)+i(2)*nz+(1+i(1))*ny*nz i(3)+(1+i(2))*nz+(1+i(1))*ny*nz 1+i(3)+(1+i(2))*nz+(1+i(1))*ny*nz];
            fields(k,:) = reshape(reshape(pb(2:4, idx, :),[24 noe])*voltages,[3 8])*coe;
        end
    end
end
