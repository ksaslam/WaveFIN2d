% EQ source Physics, HW 5 example
% MATLAB code to demonstrate use of the provided functions for solving 2d anti-plane wave propagation

% define grid

dx = 0.0025;
x1 = -1:dx:0;
x2 = 0:dx:1;

%[y1m x1m] = meshgrid(x1, x1);
%[y2m x2m] = meshgrid(x1, x2);
[y1m x1m] = meshgrid(x2, x1);
[y2m x2m] = meshgrid(x2, x2);

% grid data are stored in 6 arrays, each of which are 401x401
% v1, s11, s12 are the velocity and stress components for block 1
% v2, s21, s22 are for block 2

% initialize fields to initial data



s11 = ones(length(x1),length(x1));
%s11(abs(x1m+0.5) > 0.05) = 0.;
   % s11(abs(y1m+0.5) > 0.05) = 0.;
v1 = zeros(length(x1), length(x1));
s12 = zeros(length(x1),length(x1));
v2 = zeros(length(x1),length(x1));
s21 = zeros(length(x1),length(x1));
s22 = zeros(length(x1),length(x1));

v1= exp(-(x1m+.5).^2/.0125- (y2m-.5).^2/.0125);
v2= exp(-(x2m+.5).^2/.0125- (y2m-.5).^2/.0125);
s11= -v1;
s21= -v2;



% I use Heun's method written in low storage form (2nd order RK method)
% given system y' = f(y)
% array holding fields is y and change in fields is dy
% then algorithm for an RK step is the following:
% for i=1:nstages
%    dy = A(i)*dy+dt*f(y)
%    y = y+B(i)*dy
% Heun's method requires 2 stages with the coefficients below
% function requires passing a data structure containing this information

rk.nstages = 2;
rk.A = [0. -1.];
rk.B = [1. 0.5];

% reflection coefficients at external boundaries

rl = 0.;
rr = 0.;
rb = 0.;
rt = 0;

% wave speed, impedance, and dissipation coefficient

c = 1.;
z = 1.;
cdiss = 1.e-4/dx;

% frictional failure criteria

tau = 0.2;

% time integration info

cfl = 0.3;
dt = dx*cfl/c;
nt = 100;

% loop over time steps
% can optionally omit final 2 arguments if you want to use no dissipation and weak imposition of boundary conditions

tic
for i=1:nt
    [v1, s11, s12, v2, s21, s22] = wave_2d_step(v1, s11, s12, v2, s21, s22, dt, rk, dx, c, z, tau, rl, rr, rb, rt, true, cdiss);
end
toc

pcolor(x1m,y1m,v1)
hold on
pcolor(x2m,y2m,v2)
axis image
shading flat
colorbar
hold off
xlabel('x poistion')
ylabel('y poistion')
title(' Velocity amplitude')
%print('-dpdf','-r400', 'part3_600timestep.pdf') 