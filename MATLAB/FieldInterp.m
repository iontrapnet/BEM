function fields = FieldInterp(voltages, points)
    global xr yr zr pb
    fx=0;
    fy=0;
    fz=0;
    for i=1:size(pb,3)
        fx = fx + pb(2,:,i)*voltages(i);
        fy = fy + pb(3,:,i)*voltages(i);
        fz = fz + pb(4,:,i)*voltages(i);
    end
    for i=1:size(points)
        [ex,ey,ez]=EqualRangeInterpolation(points(i,:),xr,yr,zr,fx,fy,fz);
        fields(i,:) = [ex;ey;ez];
    end
end