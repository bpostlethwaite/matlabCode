function A = ctranspose(A);
%
%  CTRANSPOSE The transpose of the multiPsfMatrix matrix is
%             needed for matrix-vector multiply in iterative
%             restoration methods.  
%

%  K. Lee 1/31/02

if A.transpose == 0
  A.transpose = 1;
else 
  A.transpose = 0;
end