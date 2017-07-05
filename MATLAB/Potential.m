function [pot,fx,fy,fz] = Potential(triangles, charges, points)
    N = size(triangles, 3);
    for i = 1:N
        [a,b] = int_green3d_tri(points, triangles(:,:,i));
        pot(:,i) = charges(i)*a;
        fx(:,i) = charges(i)*b(:,1);
        fy(:,i) = charges(i)*b(:,2);
        fz(:,i) = charges(i)*b(:,3);
    end
    pot = sum(pot')/4/pi;
    fx = -sum(fx')/4/pi;
    fy = -sum(fy')/4/pi;
    fz = -sum(fz')/4/pi;
end
