function r = potential(triangles, charges, points)
    triangles = permute(triangles, [2 3 1]);
    N = size(triangles, 3);
    r = zeros(size(points,1),N);
    for i = 1:N
        r(:,i) = charges(i)*int_green3d_tri(points, triangles(:,:,i));
    end
    r = sum(r')/4/pi;
end