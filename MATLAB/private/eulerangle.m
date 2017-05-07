function [phi,theta,psi] = eulerangle(P0,Pz,Px)

% EULERANGLE calculates Euler angle that correspond to a local rotated
% cartesian system
%
% USE:
% [phi,theta,psi] = eulerangle(P0,Pz,Px)
%
% INPUT:
% 'P0': system origin
% 'Pz': point on z-axis
% 'Px': point on x-axis
%
% OUTPUT:
% 'phi': rotation along z-axis
% 'theta': rotation alon x-axis already rotated by phi
% 'psi': rotation alon y-axis already rotated by phi and theta
%
% NOTE:
% 'P0', 'Pz', 'Px' can be vector of the same length
% Notation according to H. Goldstein, C. Poole, J. Safko "Classical
% Mechanics", Addison Wesley, 3rd Edition, pp. 153
%
% VERSION:
% Date: 03.03.2010
% Copyright(C) 2010: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 22.03.2010: bug fix on cosarg with vector inputs

% check inputs
if (size(P0,1) ~= size(Pz,1)) || (size(P0,1) ~= size(Px,1))
    error('''P0'', ''Pz'', ''Px'' must have the same dimensions');
end

% traslation of local system
Uz = Pz-P0;
Ux = Px-P0;

% unit vectors
Uz = versor(Uz);
Ux = versor(Ux);

% angle phi from direction of cross product
phi = atan2(Uz(:,1),-Uz(:,2));

% angle theta from dot product of Uz and versor of z-axis
eta = -Uz(:,1).*sin(phi)+Uz(:,2).*cos(phi);
theta = acos(Uz(:,3));
ineg = find(eta > 0);
theta(ineg) = -theta(ineg);

% angle psi from dot product of Ux and versor of rotated x-axis
Uxr = [cos(phi) sin(phi) zeros(size(phi))];
eta = -Ux(:,1).*cos(theta).*sin(phi)+Ux(:,2).*cos(theta).*cos(phi)+Ux(:,3).*sin(theta);
cosarg = dot(Ux,Uxr,2);

% set cosarg between -1 and 1
cosarg(cosarg > 1) = 1;
cosarg(cosarg < -1) = -1;

psi = acos(cosarg);
ineg = find(eta < 0);
psi(ineg) = -psi(ineg);