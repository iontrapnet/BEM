clear;
clc
path='..\\Model\\4rod\\167634622912717531';
xr=[-0.01 0.01 10];
yr=[-0.01 0.01 10];
zr=[2.05 2.15 10];
FieldInit(path, xr, yr, zr);
%%%%%%%%%%%%%%%%%%%%%%%% Definition %%%%%%%%%%%%%%%%%%%%%%%%%
e = 1.6021766 * 10^-19;
m = 2.87363 * 10^-25;                    %Yb+离子质量
n_ions = 2;                              %ions number
Q_ion = 1 * e;                           %Yb+(ODEfun中也有)
lam = 369.5 * 10^-9;                     %cooling beam
kB = 1.38065 * 10^-23;                   %Boltzmann constant
T_total = 0.03;                           %time of simulation
d0 = 500;
d1 = 1000;
r0 = 1e-6*(d1 / sqrt(2) - d0 / 2);
v0 = 1e-6*(sqrt(2 * 500 * kB / 3 / m));
z0 = 1e-4;

%%%%%%%%%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%%%%%%%%
% eps = repmat(1e-2, 6 * n_ions, 1);
% options = odeset('RelTol',1e-2,'AbsTol',eps);
x0 = (2 * rand(n_ions, 1) - 1) * r0 / 200;
y0 = (2 * rand(n_ions, 1) - 1) * r0 / 200;
z0 = (2 * rand(n_ions, 1) - 1) * z0 / 100;
vx0 = (2 * rand(n_ions, 1) - 1) * v0 / 100;
vy0 = (2 * rand(n_ions, 1) - 1) * v0 / 200;
vz0 = (2 * rand(n_ions, 1) - 1) * v0 / 100;
ini = [x0; y0; z0; vx0; vy0; vz0]
%time = (1:T_total);
%[T,pv] = ode23(@ODEfun,time,ini,options);

%%%%%%%%%%%%%%%%%%%%%%%%  modified Euler method %%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01;
n_T = T_total / dt;
pos_vel = zeros(6 * n_ions, n_T + 1);
pos_vel(:, 1) = ini;
kt = zeros(6 * n_ions, 1);
for count = 1:n_T
    kt = Eulerfun((count - 1) * dt, pos_vel(:, count));
    pos_vel(:, count + 1) = pos_vel(:, count) + dt * kt;
end
    
for inter = n_ions:-1:1
    figure(inter);
    subplot(3, 1, 1); plot(pos_vel(inter, :)); title({['Ion ',num2str(inter)];'x'})
    subplot(3, 1, 2); plot(pos_vel(inter + n_ions, :)); title('y')
    subplot(3, 1, 3); plot(pos_vel(inter + 2 * n_ions, :)); title('z') 
end


