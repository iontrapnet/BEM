function fields = Field(voltages, points)
    global triangles cb xr yr zr pb
    if pb == 0
        [~, ex, ey, ez] = Potential(triangles, cb*voltages', points);
        fields = [ex;ey;ez]';
    else
        fx=0;
        fy=0;
        fz=0;
        for i=1:size(pb,3)
            fx = fx + pb(2,:,i)*voltages(i);
            fy = fy + pb(3,:,i)*voltages(i);
            fz = fz + pb(4,:,i)*voltages(i);
        end
        nop = size(points, 1);
        fields = zeros(nop, 3);
        for i=1:nop
            [ex,ey,ez]=EqualRangeInterpolation(points(i,:),xr,yr,zr,fx,fy,fz);
            fields(i,:) = [ex;ey;ez];
        end
    end
end