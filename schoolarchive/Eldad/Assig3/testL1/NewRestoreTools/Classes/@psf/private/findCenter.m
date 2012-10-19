function center = findCenter( PSF )
%
%       center = findCenter( PSF );
%
%  Given (possibly) several PSF images, this function sets the center 
%  of the PSFs (location of point source) to be the location of max 
%  entry of the PSFs.
%
%  Input:
%         PSF  -  cell array containing the PSF image(s)
%
%  Output:
%         center - array [row_index, col_index] containing the
%                  location of the center of the PSF.
%

%  J.Nagy 1/8/02

center = cell(size(PSF));

%
%  Note that this should work for 2D as well as 3D images.
%
for k = 1:size(PSF, 3)
  for i = 1:size(PSF, 1)
    for j = 1:size(PSF, 2)
       P = PSF{i,j,k};
       [ci, cj, ck] = ind2sub(size(P), find( P == max( P(:) ) ));
       c = [min(ci), min(cj), min(ck)];
       center{i,j,k} = c(1:length(size(P)));
    end
  end
end
