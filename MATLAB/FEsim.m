function [FEx, FEy, FEz] = FEsim(x, y, z, t, q)
wrf = 2 * pi *12;
Udc = 30; Urf = 500;
x=x.*1e2;
y=y.*1e2;
z=2.1+z.*1e2;
%dx=0.000001;dy=0.000001;dz=0.000001;
%U_harmonic = Udc/2 * (a1*X^2 + b1*Y^2 + c1*Z^2) +  Urf/2 * cos(wrf * t) * (a2*X^2 + b2*Y^2 + c2*Z^2);
Vf=[Udc,Urf.*cos(wrf.*t),0,Urf.*cos(wrf.*t),0,Udc];
FE=q*Field(Vf,[x y z]);
FEx=FE(:,1);
FEy=FE(:,2);
FEz=FE(:,3);
end