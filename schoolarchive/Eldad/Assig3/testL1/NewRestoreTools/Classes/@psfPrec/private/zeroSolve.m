function y = zeroSolve(P, b)
%
%  Overload matrix solves for a psfPrec object using
%  zero boundary conditions.
%

%  J. Nagy & K. Lee  3/7/02

PSFs = P.matdata;

nregions = size(PSFs);
imsize=size(b);
if length(imsize) == 1
  imsize = [imsize, 1, 1];
elseif length(imsize) == 2
  imsize = [imsize, 1];
end
if length(nregions) == 1
  nregions = [nregions,1,1];
elseif length(nregions) == 2
  nregions = [nregions,1];
end

rsize = ceil(imsize ./ nregions);

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

y = zeros(size(b));

if (P.transpose)

  for k = 1:size(PSFs,3)
    for i = 1:size(PSFs,1)
      for j = 1:size(PSFs,2)
        bsub = b(RIidx(i,1):RIidx(i,2), RJidx(j,1):RJidx(j,2), RKidx(k,1):RKidx(k,2));
        %M=PSFs;
        M = PSFs{i,j,k};
        y(RIidx(i,1):RIidx(i,2), RJidx(j,1):RJidx(j,2), RKidx(k,1):RKidx(k,2)) = ...
               real( ifftn( conj(M) .* fftn(bsub) ) );
       end
    end
  end

else

  for k = 1:size(PSFs,3)
    for i = 1:size(PSFs,1)
      for j = 1:size(PSFs,2)
        bsub = b(RIidx(i,1):RIidx(i,2), RJidx(j,1):RJidx(j,2), RKidx(k,1):RKidx(k,2));
        %M=PSFs;
        M = PSFs{i,j,k};
        y(RIidx(i,1):RIidx(i,2), RJidx(j,1):RJidx(j,2), RKidx(k,1):RKidx(k,2)) = ...
               real( ifftn( M .* fftn(bsub) ) );
       end
    end
  end

end

y = y(1:imsize(1), 1:imsize(2), 1:imsize(3));



