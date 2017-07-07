function [Fix, Fiy, Fiz] = Finter(x, y, z, q)
[n, ~] = size(x);
ke = 8.987551787e9;
Fix = zeros(n, 1);
Fiy = zeros(n, 1);
Fiz = zeros(n, 1);
for inter1 = 1:n
    for inter2 = 1:n
       if inter2 ~= inter1
           d = sqrt((x(inter1)-x(inter2))^2+(y(inter1)-y(inter2))^2+(z(inter1)-z(inter2))^2);
           Fix(inter1) = Fix(inter1) + ke * q^2 * (x(inter1) - x(inter2)) / d^3;
           Fiy(inter1) = Fiy(inter1) + ke * q^2 * (y(inter1) - y(inter2)) / d^3;
           Fiz(inter1) = Fiz(inter1) + ke * q^2 * (z(inter1) - z(inter2)) / d^3;
       end
    end
end
end

