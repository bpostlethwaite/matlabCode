% Question 6
clear all
close all
clc

A = [ 1,  0,  1
      2,  3,  5
      5,  3, -2
      3,  5,  4
     -1,  6,  3];
 
 B = [A, sum(A,2)];
 b = [4,-2,5,-2,1]';
 
 gamma = [0,10.^-[0,3,6,12]];
 I = eye(length(B(1,:)),length(B(1,:)));
            
 for g = gamma
 
     fprintf('\n')
     x = (B'*B + g*I)\B'*b;
     L2 = sqrt((B*x - b)'*(B*x - b));
     sol = sqrt(x'*x);
     fprintf('for gamma = %1.1e\n', g);
     fprintf('L2 norm  = %1.4f\n', L2)
     fprintf('Solution = %1.4f', sol)
     fprintf('\n')
 end