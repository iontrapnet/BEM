function dpv = Eulerfun(t, pos_vc, pb, cb)
e = 1.6021766e-19;
m = 2.87363e-25;                    %Yb+¿Î◊”÷ ¡ø                 
[nn, ~] = size(pos_vc);
n = nn / 6;
[FEx, FEy, FEz] = FEsim(pos_vc(1:n), pos_vc(n + 1:2 * n), pos_vc(2 * n + 1:3 * n), t, e, pb, cb);
[Fix, Fiy, Fiz] = Finter(pos_vc(1:n), pos_vc(n + 1:2 * n), pos_vc(2 * n + 1:3 * n),e);
%[Fcx, Fcy, Fcz] = Fcooling(pos_vc(3 * n + 1:4 * n), pos_vc(4 * n + 1:5 * n), pos_vc(5 * n + 1:6 * n));
dpv = zeros(nn,1);
for inter = 1:n
    dpv(inter) = pos_vc(3 * n + inter);
    dpv(n + inter) = pos_vc(4 * n + inter);
    dpv(2 * n + inter) = pos_vc(5 * n + inter);
    dpv(3 * n + inter) = 1 / m * (FEx(inter) + Fix(inter));%+ Fcx(inter));
    dpv(4 * n + inter) = 1 / m * (FEy(inter) + Fiy(inter)); %+ Fcy(inter));
    dpv(5 * n + inter) = 1 / m * (FEz(inter) + Fiz(inter)); %+ Fcz(inter));
end 
end
