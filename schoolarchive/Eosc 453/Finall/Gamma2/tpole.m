function [Rt] = tpole(T)

L = size(T);
P = T(L(1),:);
Tc = 125;
x = [1:L(2)]*10;
XI = 1:x(L(2));
% interpolate P here
PI = INTERP1(x,P,XI,'spline'); 
Tci = find(PI > Tc,1,'last');

Rt = (L(2)*10 - Tci)/(L(2)*10);
