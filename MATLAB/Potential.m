function [pot,fx,fy,fz] = Potential(triangles, charges, points)
    N = size(triangles, 3);
    nop = size(points, 1);
    pot = zeros(nop, 1);
    fx = zeros(nop, 1);
    fy = zeros(nop, 1);
    fz = zeros(nop, 1);
    for i = 1:N
        [a,b] = int_green3d_tri(points, triangles(:,:,i));
        pot = pot + charges(i)*a;
        fx = fx + charges(i)*b(:,1);
        fy = fy + charges(i)*b(:,2);
        fz = fz + charges(i)*b(:,3);
    end
    pot = pot'/4/pi;
    fx = -fx'/4/pi;
    fy = -fy'/4/pi;
    fz = -fz'/4/pi;
end
