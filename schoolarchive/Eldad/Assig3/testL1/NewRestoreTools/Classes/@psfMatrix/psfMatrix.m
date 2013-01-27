function A = psfMatrix(varargin)
%
%  CONSTRUCTOR FOR psfMatrix OBJECT
%
%  The psfMatrix class is based on a structure with five fields:
%    psf      - psf object 
%    matdata  - matrix data needed to do matrix vector multiplications
%    type     - character string indicating:
%                 'invariant', 'variant', 'separable'
%    boundary - character string array indicating type of boundary conditions
%               to be used.  Choices are:
%                 'zero', 'periodic', 'reflexive' (or 'neumann')
%               The default is zero.
%    transpose- indicates if the matrix has been transposed.
%
%  Calling Syntax:
%       A = psfMatrix             (returns object with empty fields)
%       A = psfMatrix(psfMatrixObj) 
%       A = psfMatrix(psfObj)
%       A = psfMatrix(psfObj, boundary)
%       A = psfMatrix(PSF)
%       A = psfMatrix(PSF, boundary)
%       A = psfMatrix(PSF, center)
%       A = psfMatrix(PSF, center, boundary)
%
%    where 
%       * psfMatrixObj is an already existing psfMatrix object
%       * psfObj       is a psf object (see help psf for more information)
%       * PSF          can be either a double array containing one PSF
%                      image, or a cell array containing one or more 
%                      PSF images
%       * boundary     is a character string indicating desired boundary
%                      condition. (see above)
%       * center       is either a double array or a cell array containing
%                      the (i,j) locations of the point sources of the 
%                      PSF images
%

%  J. Nagy & K. Lee  1/12/02


switch nargin

case 0
  A.psf = psf;
  A.matdata = [];
  A.type = '';
  A.boundary = '';
  A.transpose = 0;
  A = class(A, 'psfMatrix');

case 1
  % if single argument of class psfMatrix, return it
  % otherwise create a psfMatrix from either:
  %    - the given psfObj, or
  %    - from a double array containing the PSF image, or
  %    - from a cell array containing one or more PSF images
  %
  if ( isa( varargin{1}, 'psfMatrix' ) )
    A = varargin{1};
  elseif ( isa( varargin{1}, 'psf' ) )
    P = varargin{1};
    P = adjustPsfSize(P);
    A.psf = P;
    A.matdata = constructMatrix( P.image, P.center );
    if prod(size( P.image )) > 1
      A.type = 'variant';
    else
      A.type = 'invariant';
    end
    A.boundary = 'zero';
    A.transpose = 0;
    A = class(A, 'psfMatrix');
  elseif ( isa(varargin{1}, 'double') | isa(varargin{1}, 'cell') )
    P = psf(varargin{1});
    P = adjustPsfSize(P);
    A.psf = P;
    A.matdata = constructMatrix( P.image, P.center );
    if prod(size( P.image )) > 1
      A.type = 'variant';
    else
      A.type = 'invariant';
    end
    A.boundary = 'zero';
    A.transpose = 0;
    A = class(A, 'psfMatrix');
  else
    error('Incorrect argument type')
  end

case 2
  % create object using specific values
  if ( isa(varargin{1}, 'psf') )
    P = varargin{1};
    P = adjustPsfSize(P);
    A.psf = P;
    A.matdata = constructMatrix( P.image, P.center );
    if prod(size( P.image )) > 1
      A.type = 'variant';
    else
      A.type = 'invariant';
    end
    A.boundary = varargin{2};
    A.transpose = 0;
  elseif ( isa(varargin{1}, 'double') | isa(varargin{1}, 'cell') )
    if ( isa(varargin{2}, 'double') | isa(varargin{2}, 'cell') )
      P = psf(varargin{1}, varargin{2});
      P = adjustPsfSize(P);
      A.psf = P;
      A.matdata = constructMatrix( P.image, P.center );
      if prod(size( P.image )) > 1
        A.type = 'variant';
      else
        A.type = 'invariant';
      end
      A.boundary = 'zero';
      A.transpose = 0;
    elseif ( isa(varargin{2}, 'char') )
      P = psf(varargin{1});
      P = adjustPsfSize(P);
      A.psf = P;
      A.matdata = constructMatrix( P.image, P.center );
      if prod(size( P.image )) > 1
        A.type = 'variant';
      else
        A.type = 'invariant';
      end
      A.boundary = varargin{2};
      A.transpose = 0;
    else
      error('Incorrect argument type')
    end
  else
    error('Incorrect argument type')
  end
  A = class(A, 'psfMatrix');

case 3
  %
  P = psf(varargin{1}, varargin{2});
  P = adjustPsfSize(P);
  A.psf = P;
  A.matdata = constructMatrix( P.image, P.center );
  if prod(size( P.image )) > 1
    A.type = 'variant';
  else
    A.type = 'invariant';
  end
  A.boundary = varargin{3};
  A.transpose = 0;
  A = class(A, 'psfMatrix');

otherwise
  error('Incorrect number of input arguments.')
end



