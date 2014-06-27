clc

% Y_ml = spharm(1,1,pi/4,pi/4)

% P = sqrt(2)/2;
% N = sqrt(3/4/pi * 0.5);
% E = exp(4/pi * 1i);
% Y = P * N * E
% P*N
%P*N*i^(m+abs(m)) matches SH = gsl_sf_legendre_sphPlm(1, 1, cos(pi/4));

% A = get_SH_descriptor('A.stl',0);
% B = get_SH_descriptor('A.stl',pi/6);
whatever = get_SH_descriptor('femur.stl',0);