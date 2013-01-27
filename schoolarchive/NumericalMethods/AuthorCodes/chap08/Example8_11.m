% Example 8.11 : finding all (real) eigenvalues for a simple example

% Matrix from Example 8.9
A = [.5, -.1, -.5, .4; -.1, .3, -.2, -.3; -.3, -.2, .6, .3; .1, -.3, .3, 1];

tol = 1.e-12;
[lambda,itn] = qreig (A,tol)