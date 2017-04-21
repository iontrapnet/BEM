function r = potential(triangles, charges, points)
N = size(triangles, 3);
r = zeros(size(points,1),N);
for i = 1:N
    r(:,i) = int_green3d_tri(points, triangles(:,:,i));
end
r = r*charges./4./pi;%ALPHA*PSI
end