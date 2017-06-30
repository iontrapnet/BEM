function [Fcx, Fcy, Fcz] = Fcooling(vx, vy, vz)
lam = 369.5 * 10^-9;                     % cooling beam
hbar = 6.626070 * 10^-34;
decay = 2 * pi * 19.7 *10^6;             %decay rate
detuning = -2 * decay;
rabi = decay / sqrt(2);                  %´ý¶¨
s = 2 * abs(rabi)^2 / decay^2;    %saturation paramater
kx = 0; ky = 0; kz = 2*pi/lam;                  %wave vector
Vel = sqrt(vx.^2 + vy.^2 + vz.^2);
[n, ~] = size(vx);
pee1 = repmat(s / 2, n, 1) ./ (1 + s + (2 * (detuning - kx * vx - ky * vy - kz * vz) / decay).^2);
pee2 = repmat(s / 2, n, 1) ./ (1 + s + (2 * (detuning + kx * vx + ky * vy + kz * vz) / decay).^2);
Fcx = - hbar * 2 * pi / lam * decay * s / 2 * abs(pee1 + pee2) .* (vx ./ Vel);
Fcy = - hbar * 2 * pi / lam * decay * s / 2 * abs(pee1 + pee2) .* (vy ./ Vel);
Fcz = - hbar * 2 * pi / lam * decay * s / 2 * abs(pee1 + pee2) .* (vz ./ Vel);
end
