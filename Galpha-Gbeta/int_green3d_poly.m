function [I1,Igrad] = int_green3d_poly(Pglo,V)

% INT_GREEN3D_POLY integrates '1/r' and 'grad(1/r)' over a plane polygon
%
% USE:
% [I1,Igrad] = int_green3d_poly(Pglo,V,R)
%
% INPUTS:
% 'Pglo': points to calculate integral (in global coordinates)
% 'V': vertices of triangle
% 'integrand': '1/r', 'grad(1/r)', (case insensitive)
%
% OUTPUTS:
% 'I1': value of integral of 1/r
% 'Igrad': value of integral of Igrad(1/r)
%
% NOTE:
% See R. Graglia, "Numerical Integration of the linear shape functions
% times the 3-d Green's function or its gradient on a plane triangle", IEEE
% Transactions on Antennas and Propagation, vol. 41, no. 10, Oct 1993,
% pp. 1448--1455
%
% VERSION:
% Date: 13.12.2012
% Copyright(C) 2012-2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 16.12.2012: removed INTEGRAND from input list
% 09.09.2014: refurbished routine

% tolerance
TOL = 1e-6;

% number of field points
nfp = size(Pglo,1);

% number of vertices
nV = size(V,1);

% local edge-to-node
e2n = [1:nV; 2:nV 1]'; %[1 2; 2 3; 3 4; 4 1];

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
w0 = Ploc(:,3);

% vertices in local coordinates
Vloc = [V(:,1)-V(1,1) V(:,2)-V(1,2) V(:,3)-V(1,3)];  % translation
Vloc = Vloc*R;                      % rotation: (R'*Vloc')' = Vloc*R;

% edge lengths
L = Vloc(e2n(:,2),1:2)-Vloc(e2n(:,1),1:2);
L = sqrt(sum(L.^2,2));

% threshold for small numbers
threshold = TOL*min(L);
w0(abs(w0) < threshold) = 0;

I1 = zeros(nfp,1);
Igradloc = zeros(nfp,3);

for i = 1:nV
    % versors normal to edges Fig. 1(b)
    m = [Vloc(e2n(i,2),2)-Vloc(e2n(i,1),2) Vloc(e2n(i,1),1)-Vloc(e2n(i,2),1) 0];
    m = versor(m);
    
    % calculate parameter p of intersection
    % see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    p = -((Vloc(e2n(i,1),1)-Ploc(:,1))*(Vloc(e2n(i,2),1)-Vloc(e2n(i,1),1))+...
        (Vloc(e2n(i,1),2)-Ploc(:,2))*(Vloc(e2n(i,2),2)-Vloc(e2n(i,1),2)))/L(i)^2;
    
    % intersection points
    Pintersect = [Vloc(e2n(i,1),1)+(Vloc(e2n(i,2),1)-Vloc(e2n(i,1),1))*p...
        Vloc(e2n(i,1),2)+(Vloc(e2n(i,2),2)-Vloc(e2n(i,1),2))*p];
    
    % useful quantities for integration
    sminus = -p*L(i);
    splus = (1-p)*L(i);
    
    tplus = sqrt((Ploc(:,1)-Vloc(e2n(i,2),1)).^2+(Ploc(:,2)-Vloc(e2n(i,2),2)).^2);
    tminus = sqrt((Ploc(:,1)-Vloc(e2n(i,1),1)).^2+(Ploc(:,2)-Vloc(e2n(i,1),2)).^2);
    
    t0vec = Ploc(:,1:2)-Pintersect;
    isgn = -sign(t0vec(:,1)*m(1)+t0vec(:,2)*m(2));
    t0 = isgn.*sqrt(sum(t0vec.^2,2));
    
    R0 = sqrt(t0.^2+w0.^2);                    % line 1 pp. 1450
    Rminus = sqrt(tminus.^2+w0.^2);            % line 2 pp. 1450
    Rplus = sqrt(tplus.^2+w0.^2);              % line 2 pp. 1450
    
    f2 = log((Rplus+splus)./(Rminus+sminus));              % eq. (11)
    beta = atan(t0.*splus./(R0.^2+abs(w0).*Rplus))-...
        atan(t0.*sminus./(R0.^2+abs(w0).*Rminus));      % eq. (14)
    
    % field point in the plane of triangle and aligned with edges
    id = (abs(w0) < threshold) & (abs(t0) < threshold);
    beta(id) = 0;
    f2(id) = abs(log(splus(id)./sminus(id)));                  % abs(lim t->0 eq. 15)
    f2(~isfinite(f2)) = 0;                                        % fix value for point on triangle corners (undocumented)
    
    I1 = I1+t0.*f2-abs(w0).*beta;         % integral value eq. (19)
    Igradloc = Igradloc - [m(1)*f2 m(2)*f2 sign(w0).*beta];       % integral value eq. (34)
    
end
Igrad = Igradloc*R';

end