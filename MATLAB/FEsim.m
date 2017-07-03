function [FEx, FEy, FEz] = FEsim(x, y, z, t, q)
global xr yr zr pb
wrf = 2 * pi *12;
Udc = 30; Urf = 500;
x=x.*1e2;
y=y.*1e2;
z=2.1+z.*1e2;
[x y z]
%dx=0.000001;dy=0.000001;dz=0.000001;
%U_harmonic = Udc/2 * (a1*X^2 + b1*Y^2 + c1*Z^2) +  Urf/2 * cos(wrf * t) * (a2*X^2 + b2*Y^2 + c2*Z^2);
Vf=[Udc,Urf.*cos(wrf.*t),0,Urf.*cos(wrf.*t),0,Udc];
fx=0;
fy=0;
fz=0;
for i=1:size(pb,3)
    fx = fx + pb(2,:,i)*Vf(i);
    fy = fy + pb(3,:,i)*Vf(i);
    fz = fz + pb(4,:,i)*Vf(i);
end
for i=1:size(x)
    [ex,ey,ez]=EqualRangeInterpolation([x(i) y(i) z(i)],xr,yr,zr,fx,fy,fz);
    FEx(i)=ex;
    FEy(i)=ey;
    FEz(i)=ez;
end
FEx=q*FEx;
FEy=q*FEy;
FEz=q*FEz;
end