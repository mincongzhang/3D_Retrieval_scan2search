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
clc
clear
[Ymn,PHI,THETA,Xm,Ym,Zm]=spharm_array(2,2);