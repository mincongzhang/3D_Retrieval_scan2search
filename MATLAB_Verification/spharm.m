function Ymn = spharm(L,M,THETA,PHI)
if nargin==0
  L=2;   % DEGREE
  M=1;   % ORDER
  THETA = pi/4;
  PHI = pi/4;
end

if L<M, error('The ORDER (M) must be less than or eqaul to the DEGREE(L).'); end

%THETA  Azimuthal/Longitude/Circumferential
%PHI    Altitude /Latitude /Elevation

Lmn=legendre(L,cos(PHI));

if L~=0
  Lmn=squeeze(Lmn(abs(M)+1,:,:));
end

a1=((2*L+1)/(4*pi));
a2=factorial(L-M)/factorial(L+M);
C=sqrt(a1*a2);

Ymn=C*Lmn.*exp(i*M*THETA);
end