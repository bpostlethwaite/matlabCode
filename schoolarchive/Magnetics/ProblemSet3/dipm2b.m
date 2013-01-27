function [bx,by,bz] = dipm2b(XX,YY,ZZ,mvec)

% Inputs:   XX, YY, ZZ:  arrays giving x,y,z position
%           of observer relative to the dipole:  can be
%           1D or 2D.  These should be in METERS.
%           mvec:   x,y,z components of the magnetic dipole
%           Should be in Am
% Outputs:  bvec=(bx,by,bz) - will be in TESLAS

% edit to input magnetization eventually?

mu0=4*pi*1e-07;     % SI units
cc=mu0/(4*pi);

r=sqrt(XX.^2+YY.^2+ZZ.^2);
r2=r.^2;
r5=r.^5;

mdotr=XX.*mvec(1)+YY.*mvec(2)+ZZ.*mvec(3);
bx=cc*(3.*mdotr.*XX-r2.*mvec(1))./r5;
by=cc*(3.*mdotr.*YY-r2.*mvec(2))./r5;
bz=cc*(3.*mdotr.*ZZ-r2.*mvec(3))./r5;

