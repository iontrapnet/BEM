function fields = FieldCalc(voltages, points)
    global triangles cb
    [~, ex, ey, ez] = Potential(triangles, cb*voltages', points);
    fields = [ex;ey;ez]';
end