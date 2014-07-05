
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
% P = sqrt(2)/2;
% N = sqrt(3/4/pi * 0.5);
% m = 0;
% l = 0;
% SH_manual = P*N*i^(m+abs(m))
% 
% %MATLAB

R_vector = [1:32];
phi_vector = ones(1,32);
theta_vector = ones(1,32)*3;
for r = 1:32
    for l = 1:32
        for m = -l:l
            if(m>=0)
                Y_ml(l,r) = spharm(l,m,phi_vector(r),theta_vector(r));
            else
                Y_temp = spharm(l,-m,phi_vector(r),theta_vector(r));
                Y_ml(l,r) = (-1)^(-m) * conj(Y_temp);
            end
        end
    end
end
figure,
bar3(abs(Y_ml))


%%
% A = get_SH_descriptor('A.stl',0);
% B = get_SH_descriptor('A.stl',pi/6);
% B = get_SH_descriptor('femur.stl',pi/6);