function P = ctranspose(P);
%
%  CTRANSPOSE The transpose of the psfPrec matrix is
%             needed for matrix-vector multiply in iterative
%             restoration methods.  
%

%  J. Nagy  7/8/01

if P.transpose == 0
  P.transpose = 1;
else 
  P.transpose = 0;
end
