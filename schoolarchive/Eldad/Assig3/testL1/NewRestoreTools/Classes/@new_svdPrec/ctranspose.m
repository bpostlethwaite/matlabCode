function P = ctranspose(P);
%
%  CTRANSPOSE The transpose of the svdPrec matrix is
%             needed for matrix-vector multiply in iterative
%             restoration methods.  
%

%  J. Nagy  6/22/02

U = P.u;
P.u = P.v;
P.v = U;
P.s = conj(P.s);