%% PROBLEM SET 3:  QUESTION 4

clear all;
close all;

t2nt=1e09;            % Teslas to nT

% Set up dipole moment, magnitude mm, declination D, inclination I 



% Compute mvec from mm, D, I



% Set up dipole position, xq, yq, zq



% take observations at positions xp,yp,zp



% need dipole-observer vector, e.g. x=xp-xq


% Use meshgrid to make a 2-D array of x,y values and also make a 2-D array
% of z values (z will of course be the same at every point in the grid).
% Call these 2D arrays XX, YY, ZZ



%  Call the function dipm2b



%  Convert the Bx, By, Bz to nT and calculate Bh



% Make fig with 4 suplots showing contour maps of Bx, By, Bz, Bh (Bh = horizontal comp mag
% field).  MAKE sure to plot your contour map in terms of your ACTUAL
% observation positions (i.e. xp, rp) NOT the relative position of the observations wrt
% the dipole (x,y). 
