% Example 9.2 : finding roots for Example 9.1 using newtons

% First initial guess
x0 = [1,1]';
%find root
[x,k] = newtons(@func,x0,1e-6,20)

% Second initial guess
x0 = [-1,1]';
%find root
[x,k] = newtons(@func,x0,1e-6,20)
