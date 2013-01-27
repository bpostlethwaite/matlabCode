function y = mldivide(arg1, arg2)
%
%  Overload backslash operation for psfPrec object.
%
%  Implements P \ b and P' \ b   for psfPrec object
%  P and vector b.  Result is returned as a vector y. 
%

%  J. Nagy  7/2/01

if ( isa(arg1, 'psfPrec') )
  %
  % See if arg2 is vec(image), and if so, reshape to be image.
  %
  if prod(size(arg2)) == length(arg2)
    rs = true;
    P = arg1.matdata;
    P = P{1,1};
    if prod(size(P)) == length(arg2)
      arg2 = reshape(arg2, size(P));
    else
      error('Cannot determine reshape size from PSF size.')
    end
  else
    rs = false;
  end

  switch arg1.boundary

  case 'zero'
    y = zeroSolve(arg1, arg2);

  case 'periodic'
    y = periodicSolve(arg1, arg2);

  case 'neumann'
    y = neumannSolve(arg1, arg2);

  otherwise
    error('Invalid boundary condition')
  end
  
  %
  %  Check to see if input arg2 was vec(image).  If so, reshape
  %  output so that it is a vec as well.
  if rs 
    y = y(:);
  end 

else

  error('incorrect argument type')

end
