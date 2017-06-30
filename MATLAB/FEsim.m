function [FEx, FEy, FEz] = FEsim(x, y, z, t, q)
global triangles charge noe
wrf = 2 * pi *12;
Udc = 30; Urf = 500;
x=x.*1e4;
y=y.*1e4;
z=z.*1e4+21;%ºÁÃ×
%dx=0.000001;dy=0.000001;dz=0.000001;
%U_harmonic = Udc/2 * (a1*X^2 + b1*Y^2 + c1*Z^2) +  Urf/2 * cos(wrf * t) * (a2*X^2 + b2*Y^2 + c2*Z^2);
chargei=0;
Vf=[Udc,Urf.*cos(wrf.*t),0,Urf.*cos(wrf.*t),0,Udc];
for i=1:noe
    chargei=chargei+charge(:,i)*Vf(i);
end
points=[x,y,z];
[~,FEx,FEy,FEz]= Potential(triangles, chargei, points);
FEx=q*FEx;
FEy=q*FEy;
FEz=q*FEz;
end