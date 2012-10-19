function y = mtimes(arg1, arg2)
%
%   Overload Matrix multiplication operations for psfMatrix
%
%   Implement A*x and A'*x for psfMatrix object A and "vector" x.  
%   (x is either a 2-d or 3-d image).
%
%   Result is returned as a vector y.  
%
%   If the blur is spatially variant, we use piecewise constant 
%   interpolation of the individual PSFs.
%

%  J. Nagy & K. Lee  1/29/02
%  
%  Modifications: 
%    J. Nagy - modify so that a vector version of a 2-d image
%              can be used for arg2.  The vector should be a
%              column ordered version of the image; that is,
%              x = vec(X) = X(:)

if ( isa( arg1, 'psfMatrix' ) )
  % check to see of arg2 is a scalar
  if (( isa ( arg2, 'double' )) & (length(arg2)==1))
    y=arg1;
    M = y.matdata;
    for k = 1:size(M,3);
      for j = 1:size(M,2);
        for i = 1:size(M,1)
           M{i,j,k} = arg2*M{i,j,k};
        end
      end
    end
    y.matdata = M;

  else
    
    %
    %  Check to see if arg2 is vec(image), and if so reshape
    %  to be image.
    %
    if prod(size(arg2)) == length(arg2)
      rs = true;
      if prod(size(arg1.psf)) == length(arg2)
        arg2 = reshape(arg2, size(arg1.psf));
      else
        error('Cannot determine reshape size from PSF size.')
      end
    else
      rs = false;
    end
    
    psfSize = size(arg1.psf);
    padsize = round(size(arg1.psf)/2);

    switch arg1.boundary
    case 'zero'
      X = padarray(arg2, padsize, 0, 'both');
    case 'reflexive'
      X = padarray(arg2, padsize, 'symmetric', 'both');
    case 'periodic'
      X = padarray(arg2, padsize, 'circular', 'both');
    case 'neumann'
      X = padarray(arg2, padsize, 'symmetric', 'both');
    case 'replicate'
      X = padarray(arg2, padsize, 'replicate', 'both');
    otherwise
      error('Illegal boundary condition');
    end  
  
    switch arg1.type
    case 'invariant'

      psfMatData = arg1.matdata;
      if (arg1.transpose)
        y = invariantMultiply(conj(psfMatData{1}), X, padsize);
      else
        y = invariantMultiply(psfMatData{1}, X, padsize);
      end

    case 'variant'

      psfMatData = arg1.matdata;
      if (arg1.transpose)
        y = variantTransposeMultiply(psfMatData, X, padsize);
      else
        y = variantMultiply(psfMatData, X, padsize);
      end

    otherwise
      error('Invalid psfMatrix type')
    end
    
    %
    %  Check to see if input arg2 was vec(image).  If so, reshape 
    %  output so that it is a vec as well.
    if rs
      y = y(:);
    end
    
  end

elseif (( isa(arg1, 'double')) & (length(arg1)==1))
  % check to see if arg1 is a scalar
  %%y=arg2;
  %%M = y.matdata;
  %%for k = 1:size(M,3);
  %%  for j = 1:size(M,2);
  %%    for i = 1:size(M,1);
  %%      M{i,j,k} = arg1*M{i,j,k};
  %%    end
  %%  end
  %%end
  %%y.matdata = M;
  
  y=arg2;
  M = y.matdata;
  A = y.psf;, A1 = A.image;
  for k = 1:size(M,3);
    for j = 1:size(M,2);
      for i = 1:size(M,1);
        M{i,j,k} = arg1*M{i,j,k};
        A1{i,j,k} = arg1*A1{i,j,k};
      end
    end
  end
  y.matdata = M;
  A.image = A1;
  y.psf = A;
 

    
 %   y=arg2;
 %   psfmatdata=arg2.matdata;
 %   y.matdata=arg1*(psfmatdata{1});
  
else
  error('Right multiplication is not implemented.')
end
