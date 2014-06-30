clc

% Y_ml = spharm(1,1,pi/4,pi/4)
% Y_ml(ml_count,1) = spharm(idx_l,idx_m,theta(idx_n),phi(idx_n));

% P = sqrt(2)/2;
% N = sqrt(3/4/pi * 0.5);
% E = exp(4/pi * 1i);
% Y = P * N * E
% P*N
%P*N*i^(m+abs(m)) matches SH = gsl_sf_legendre_sphPlm(1, 1, cos(pi/4));

%%
%TEST SH in C++
%C++: 
%SH = gsl_sf_legendre_sphPlm(1, 1, cos(pi/4));

%manual
P = sqrt(2)/2;
N = sqrt(3/4/pi * 0.5);
m = 1;
l = 1;
SH_manual = P*N*i^(m+abs(m))

%MATLAB
SH_matlab = spharm(l,m,pi/4,pi/4)/exp(1i*m*pi/4)



%%
% A = get_SH_descriptor('A.stl',0);
% B = get_SH_descriptor('A.stl',pi/6);
% whatever = get_SH_descriptor('femur.stl',0);