function y=mtimes(arg1,arg2)
%
%  Overload matrix multiplication operations for multiPsfMatrix
%
%  Implement A*x, x'*A, A'*x, x'*A' for multiPsfMatrix object
%  A and vector x.  Result is returned as a vector y.  Here we 
%  use piecewise constant interpolation if the psfMatrix is 
%  spatially variant.
%

%  K. Lee 1/31/02

if (isa ( arg1, 'multiPsfMatrix'))
  K=arg1.psfMatrices;
  L=length(K);
  if ( arg1.transpose == 1)
    [mm,n] = size(arg2);
    m = mm/L;
    if ( fix(m) ~= m )
      error('Incorrect input arguments')
    end
    y = zeros(m,n);
    for i = 1:L
      B=K{i};
      %B=B';
      y = y + B' * arg2((i-1)*m+1:i*m,:);
    end
  else
    [m,n] = size(arg2);
    y = zeros(m*L,n);
    for i = 1:L
      y((i-1)*m+1:i*m,:) = K{i} *arg2;
    end
  end
else
  error('Right multiplication is not yet implemented')
end