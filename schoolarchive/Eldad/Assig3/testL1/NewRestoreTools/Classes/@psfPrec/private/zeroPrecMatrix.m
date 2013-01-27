function psfPrecMatData = zeroPrecMatrix( PSFs, center, b, tol )
%
%       psfPrecMatData = zeroPrecMatrix( PSFs, center, b, tol );
%
%  Construct psfPrec matrix data using zero boundary condtions.
%
%  Given (several) PSFs, the locations of the corresponding point sources, 
%  and the boudary condition, this function sets up the data (eigenvalues)
%  needed for a circulant preconditioner defined by the PSFs.
%
%  Input:
%        PSFs  -  cell array containing the PSFs
%      center  -  cell array with {i,j,k} entry containing
%                 [row_index, col_index, k_index] of point source location for
%                 PSF{i,j,k}
%           b  -  blurred image
%
%  Optional Input:
%         tol  -  cell array containing the tolerance used to "regularize" the 
%                 preconditioner for each individual PSF
%                 (e.g., using Hanke, Nagy, Plemmons approach) 
%
%  Output:
%   psfPrecMatData -  cell array containing the data needed for the
%                     preconditioner solve routines
%

%  J. Nagy  &  K. Lee 2/8/02
if nargin < 4
  tol = [];
end

imsize=size(b);
if length(imsize)==1
  imsize = [imsize,1,1];
elseif length(imsize) == 2
  imsize = [imsize,1];
end
nregions=size(PSFs);
if length(nregions)==1
  nregions = [nregions,1,1];
elseif length(nregions) == 2
  nregions = [nregions,1];
end
rsize = ceil(imsize ./ nregions);

%
%  In order for this to be consistent for 2-D and 3-D images, we need to make
%  sure there is a third dimension ...
%
if length(imsize) == 1
  imsize = [imsize, 1, 1];
  rsize = [rsize, 1, 1];
  nregions = [nregions, 1, 1];
elseif length(imsize) == 2
  imsize = [imsize, 1];
  rsize = [rsize, 1];
  nregions = [nregions, 1];
end

%
%  Coding the rest of this will be easier if all of the image subregions
%  have the same dimensions.  If it's not, we pad with a few zeros to make
%  it so ...
%
padsize1 = rsize .* nregions - imsize;
if any( padsize1 < 0 )
  error('Something is wrong here ...')
end
b = padarray(b, padsize1, 'post');

%
%  Now we get information about beginning and ending indices of subregions
%  so we can "get" subregions correctly ...
%
[RIidx, RJidx, RKidx] = region_indices( nregions, rsize );

psfPrecMatData = cell(nregions);

if isempty(tol)

  for k = 1:size(PSFs,3)
    for i = 1:size(PSFs,1)
      for j = 1:size(PSFs,2)
        PSF = PSFs{i,j,k};
        c = center{i,j,k};
        bsub = b(RIidx(i,1):RIidx(i,2), RJidx(j,1):RJidx(j,2), RKidx(k,1):RKidx(k,2));
        psfPrecMatData{i,j,k} = zeroPrecMatrixOnePsf(PSF, c, bsub);
      end
    end
  end

elseif isnumeric(tol)                            % if a single value for tol is given, 
   % tol = num2cell( tol * ones(size(PSFs)) );  % use it for all PSFs
  

  for k = 1:size(PSFs,3)
    for i = 1:size(PSFs,1)
      for j = 1:size(PSFs,2)
        PSF = PSFs{i,j,k};
        c = center{i,j,k};
        %t = tol{i,j,k};
        bsub = b(RIidx(i,1):RIidx(i,2), RJidx(j,1):RJidx(j,2), RKidx(k,1):RKidx(k,2));
        psfPrecMatData{i,j,k} = zeroPrecMatrixOnePsf(PSF, c, bsub, tol);
      end
    end
  end

else

  for k = 1:size(PSFs,3)
    for i = 1:size(PSFs,1)
      for j = 1:size(PSFs,2)
        PSF = PSFs{i,j,k};
        c = center{i,j,k};
%        t = 'help';
        bsub = b(RIidx(i,1):RIidx(i,2), RJidx(j,1):RJidx(j,2), RKidx(k,1):RKidx(k,2));
        psfPrecMatData{i,j,k} = zeroPrecMatrixOnePsf(PSF, c, bsub, tol);
      end
    end
  end
end



