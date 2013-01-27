function center = findCenter( PSF )
%
%       center = findCenter( PSF );
%
%  Given a single PSF image, this function sets the center of the PSF 
%  (location of point source) to be the max  entry of the PSF.
%
%  Input:
%         PSF  -  double array containing the PSF image
%
%  Output:
%         center - array [row_index, col_index] containing the
%                  location of the center of the PSF.
%

%  J. Nagy  2/25/01

[row_index, col_index] = find( PSF == max(max(PSF)) );
center = [row_index, col_index];