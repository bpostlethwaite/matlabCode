% Example 9.13 : a small linear programming problem

A = [2,1,1,0,0; 1,1,0,1,0; 1,0,0,0,1];
b = [75;60;25];
c = [-150,-100,0,0,0]';

[x,gap,nbas] = lpm (A,b,c)