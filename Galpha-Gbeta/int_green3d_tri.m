function [I1,Igrad] = int_green3d_tri(Pglo,V)

% INT_GREEN3D_TRI integrates '1/r' and 'grad(1/r)' over a triangle
%
% USE:
% I = int_green3d_TRI(Pglo,V,integrand)
%
% INPUTS:
% 'Pglo': point(s) to calculate integral (in global 3d coordinates)
% 'V': vertices of triangle (in 3d)
%
% OUTPUTS:
% 'I1': value of the integral of '1/r'
% 'Igrad': value of the integral of 'grad(1/r)'
%
% NOTE:
% See R. Graglia, "Numerical Integration of the linear shape functions
% times the 3-d Green's function or its gradient on a plane triangle", IEEE
% Transactions on Antennas and Propagation, vol. 41, no. 10, Oct 1993,
% pp. 1448--1455
%
% VERSION:
% Date: 03.03.2010
% Copyright(C) 2010-2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 22.03.2010: introduced threshold for check w0,t0 == 0
% 29.01.2010: Changed local orientation. Basically:
% old -> new
%  1  ->  3
%  2  -> -2
%  3  ->  1
% 23.05.2011: bug fix on dimension of F2 and BETA matrices
% 16.01.2012: bug fix when checking w0 == 0
% 06.07.2012: bug fix for points on triangle corners
% 10.09.2012: moved W0REP calculation after threshold application
% 09.09.2014: removed INTEGRAND input
% 09.09.2014: refurbished routine

% number of field points
nfp = size(Pglo,1);

% local-to-global rotation matrix
P0 = V(1,:);         % origin of local coordinate system
Px = V(2,:);         % point on local x-axis
v1 = V(2,:)-V(1,:);
v2 = V(3,:)-V(1,:);
dz = versor([v1(2)*v2(3)-v2(2)*v1(3) v2(1)*v1(3)-v1(1)*v2(3) v1(1)*v2(2)-v2(1)*v1(2)]);   % dz = cross(v1,v2,2), versor to avoid small areas
Pz = P0+dz;                                                                               % point on local z-axis: Pz = P0+dz
[phi,theta,psi] = eulerangle(P0,Pz,Px);                                                   % Euler's angles
R = rotation(phi,theta,psi);

% field points in local coordinates
Ploc = Pglo-ones(nfp,1)*V(1,:); % repmat2(V(1,:),nfp,1);          % translation of field points Pglo
Ploc = Ploc*R;                              % rotation: (R'*Ploc')' = Ploc*R
u0 = Ploc(:,1);                             % notation according to Graglia's paper
v0 = Ploc(:,2);
w0 = Ploc(:,3);
%w0rep = w0(:,ones(3,1));

% vertices in local coordinates
Vloc = V-[V(1,:); V(1,:); V(1,:)];  % translation
Vloc = Vloc*R;                      % rotation: (R'*Vloc')' = Vloc*R;
l3 = Vloc(2,1);                     % notation according to Graglia's paper
u3 = Vloc(3,1);
v3 = Vloc(3,2);

% edge lengths
l1 = sqrt((l3-u3)^2+v3^2);
l2 = sqrt(u3^2+v3^2);

% threshold for small numbers
threshold = 1e-6*min([l1,l2,l3]);
w0(abs(w0) < threshold) = 0;
w0rep = w0(:,ones(3,1));

% versors normal to edges Fig. 1(b)
m = [              % m = s x w, with w = [0 0 1]
    v3 l3-u3 0     % m1
    -v3 u3 0       % m2
    0 -l3 0        % m3
    ];
m = versor(m);

% useful quantities for integration
sminus = zeros(nfp,3);                        % eq. (3)
sminus(:,1) = -((l3-u3)*(l3-u0)+v3*v0)/l1;
sminus(:,2) = -(u3*(u3-u0)+v3*(v3-v0))/l2;
sminus(:,3) = -u0;

splus = zeros(nfp,3);                         % eq. (3)
splus(:,1) = ((u3-l3)*(u3-u0)+v3*(v3-v0))/l1;
splus(:,2) = (u3*u0+v3*v0)/l2;
splus(:,3) = l3-u0;

t0 = zeros(nfp,3);                            % eq. (4)
t0(:,1) = ((u3-l3)*v0+v3*(l3-u0))/l1;
t0(:,2) = (v3*u0-u3*v0)/l2;
t0(:,3) = v0;

tplus = zeros(nfp,3);                         % eq. (5)
tplus(:,1) = sqrt((u3-u0).^2+(v3-v0).^2);
tplus(:,2) = sqrt(u0.^2+v0.^2);
tplus(:,3) = sqrt((l3-u0).^2+v0.^2);

tminus = zeros(nfp,3);                        % eq. (5)
tminus(:,1) = tplus(:,3);
tminus(:,2) = tplus(:,1);
tminus(:,3) = tplus(:,2);

R0 = sqrt(t0.^2+w0rep.^2);                    % line 1 pp. 1450

Rminus = sqrt(tminus.^2+w0rep.^2);            % line 2 pp. 1450

Rplus = sqrt(tplus.^2+w0rep.^2);              % line 2 pp. 1450

% field point not in the plane of triangle
[ir,jc] = find(abs(w0rep)  >= threshold);
id1 = sub2ind(size(w0rep),ir,jc);

% field point in the plane of triangle but not aligned with edges
[ir,jc] = find((abs(w0rep) < threshold) & (abs(t0) >= threshold));
id2 = sub2ind(size(w0rep),ir,jc);

% field point in the plane of triangle and aligned with edges
[ir,jc] = find((abs(w0rep) < threshold) & (abs(t0) < threshold));
id3 = sub2ind(size(w0rep),ir,jc);

% initialize f2 and beta matrices
f2 = zeros(size(w0rep));
beta = zeros(size(w0rep));

f2(id1) = log((Rplus(id1)+splus(id1))./(Rminus(id1)+sminus(id1)));              % eq. (11)
beta(id1) = atan(t0(id1).*splus(id1)./(R0(id1).^2+abs(w0rep(id1)).*Rplus(id1)))-...
    atan(t0(id1).*sminus(id1)./(R0(id1).^2+abs(w0rep(id1)).*Rminus(id1)));      % eq. (14)

f2(id2) = log((tplus(id2)+splus(id2))./(tminus(id2)+sminus(id2)));                  % eq. (15)
beta(id2) = atan(splus(id2)./t0(id2))-atan(sminus(id2)./t0(id2));                   % eq. (17)

beta(id3) = 0;
f2(id3) = abs(log(splus(id3)./sminus(id3)));                  % abs(lim t->0 eq. 15)
f2(~isfinite(f2)) = 0;                                        % fix value for point on triangle corners (undocumented)

% integral value of '1/r'
I1 = sum(t0.*f2-abs(w0rep).*beta,2);         % integral value eq. (19)
% integral value of grad(1/r)
Igradloc = [-m(1,1)*f2(:,1)-m(2,1)*f2(:,2)-m(3,1)*f2(:,3) -m(1,2)*f2(:,1)-m(2,2)*f2(:,2)-m(3,2)*f2(:,3) -sign(w0).*sum(beta,2)];       % integral value eq. (34)
Igrad = Igradloc*R';

end


