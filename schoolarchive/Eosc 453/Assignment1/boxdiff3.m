function  dmass  = boxdiff3(t,mass,t_yr,e_GtC_yr,modt,mod)
% boxdiff - sets up a system of equations for usage by an ODE solver
global K    %bring in K from main script

%%Creates e CO2 forcing stops at 2100: note that need some zeros close in time to
%%approximate stepwise shutoff with linear interp
e = interp1(t_yr , e_GtC_yr , t , 'linear');
m = interp1(modt,mod,t);

%Create system of equations
dmass1 = K(2,1)*mass(2)+K(5,1)*mass(5)+K(6,1)*mass(6)+K(9,1)*mass(9)+K(8,1)*mass(8)...
    -K(1,2)*mass(1)-K(1,5)*mass(1)+K(7,1)*mass(7)+m*e;
dmass2 = K(1,2)*mass(1)-K(2,1)*mass(2)+K(4,2)*mass(4)+K(7,2)*mass(7)-K(2,4)*mass(2)...
    -K(2,3)*mass(2)+K(3,2)*mass(3);
dmass3 = K(2,3)*mass(2)-K(3,2)*mass(3)-K(3,4)*mass(3);
dmass4 = K(3,4)*mass(3)-K(4,2)*mass(4)+K(2,4)*mass(2);
dmass5 = K(1,5)*mass(1)-K(5,1)*mass(5)-K(5,6)*mass(5)-K(5,7)*mass(5);
dmass6 = K(5,6)*mass(5)-K(6,1)*mass(6)-K(6,7)*mass(6);
dmass7 = K(5,7)*mass(5)+K(6,7)*mass(6)-K(7,2)*mass(7)-K(7,1)*mass(7)-K(7,8)*mass(7)-K(7,9)*mass(7);
dmass8 = K(7,8)*mass(7)-K(8,1)*mass(8);
dmass9 = K(7,9)*mass(7)-K(9,1)*mass(9);

%set dmass to a column vector consisting of equations
dmass=[dmass1;dmass2;dmass3;dmass4;dmass5;dmass6;dmass7;dmass8;dmass9];


end

