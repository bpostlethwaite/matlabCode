function psfPrecMatData = constructMatrix( PSFs, center, boundary, b, tol )
%
%       psfPrecMatData = constructMatrix( PSFs, center, boundary, b, tol );
%
%  Given (several) PSFs, the locations of the corresponding point sources, 
%  and the boudary condition, this function sets up the data (eigenvalues)
%  needed for a circulant preconditioner defined by the PSFs.
%
%  Input:
%         PSF  -  cell array containing the PSFs
%      center  -  cell array with {i,j} entry containing
%                 [row_index, col_index] of point source location for
%                 PSF{i,j}
%    boundary  -  character string indicating which boundary condition
%                 is to be implemented
%           b  -  blurred image
%
%  Optional Input:
%         tol  -  tolerance used to "regularize" the preconditioner
%                 (e.g., using Hanke, Nagy, Plemmons approach) 
%
%  Output:
%   psfPrecMatData -  cell array containing the data needed for the
%                     preconditioner solve routines
%

%  J. Nagy  12/31/01

if nargin == 4
  tol = [];
end

switch boundary

case 'zero'

  psfPrecMatData = zeroPrecMatrix( PSFs, center, b, tol );

case 'periodic'

  psfPrecMatData = periodicPrecMatrix( PSFs, center, b, tol );

case 'neumann'

  psfPrecMatData = neumannPrecMatrix( PSFs, center, b, tol );

otherwise

  error('Incorrect boundary condition')

end

