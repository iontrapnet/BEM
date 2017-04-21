function R = rotation(phi,theta,psi)

% ROTATION build rotation matrix Pglo = Ploc*R' (or Ploc = Pglo*R)
% 
% USE:
% R = rotation(phi,theta,psi)
% 
% INPUT:
% 'phi': rotation along z-axis
% 'theta': rotation alon x-axis already rotated by phi
% 'psi': rotation alon y-axis already rotated by phi and theta
% 
% OUTPUT:
% 'R' rotation matrix
%
% NOTE:
% 'phi', 'theta', 'psi' can be vector of the same length
% Notation according to H. Goldstein, C. Poole, J. Safko "Classical
% Mechanics", Addison Wesley, 3rd Edition, pp. 153
%
% VERSION:
% Date: 03.03.2010
% Copyright(C) 2010: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
%

% check inputs
if (length(phi) ~= length(theta)) || (length(phi) ~= length(psi))
    error('''phi'', ''theta'', ''psi'' must have the same length');
end

cosphi = cos(phi);
sinphi = sin(phi);
costheta = cos(theta);
sintheta = sin(theta);
cospsi = cos(psi);
sinpsi = sin(psi);

R = zeros(3,3,length(phi));
R(1,1,:) = cospsi.*cosphi-costheta.*sinphi.*sinpsi;
R(1,2,:) = -sinpsi.*cosphi-costheta.*sinphi.*cospsi;
R(1,3,:) = sintheta.*sinphi;

R(2,1,:) = cospsi.*sinphi+costheta.*cosphi.*sinpsi;
R(2,2,:) = -sinpsi.*sinphi+costheta.*cosphi.*cospsi;
R(2,3,:) = -sintheta.*cosphi;

R(3,1,:) = sinpsi.*sintheta;
R(3,2,:) = cospsi.*sintheta;
R(3,3,:) = costheta;

