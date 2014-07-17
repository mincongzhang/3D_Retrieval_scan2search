% clc
% clear
% close all
% r = 1;
% idx = 1;
% theta = pi/4;
% phi = pi/2;
% for l = 5:5
%     temp = 0;
%     for m = -l:l
%         if(m>=0)
%             SH = spharm(l,m,theta,phi)
%             temp = temp + SH;
%             %Y_ml(l,r) =Y_ml(l,r) + spharm(l,m,pi/4,pi/4);
%         else
%             SH = (-1)^(-m) * conj(spharm(l,-m,theta,phi))
%             temp = temp + SH;
%             %Y_temp = spharm(l,-m,pi/4,pi/4);
%             %Y_ml(l,r) =Y_ml(l,r) + (-1)^(-m) * conj(Y_temp);
%         end
% 
%     end
%     Y_ml = temp;
%     finalX = real(temp);
%     finalY = imag(temp);
% end
l = 1;
m = 0;
theta = pi/4;
phi = pi/4;
l = 2;
YML_star = 0;
for m = -l:l
    YML_star = conj(spharm(l,m,theta,phi));
end

RECOVER = 0;
for t = 0:0.01:pi
    for p = 0:0.01:2*pi
        RECOVER = RECOVER + YML_star
    end
end
YML_recover = 