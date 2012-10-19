function precMatData = zeroPrecOnePsf(PSF, center, b, tol)
%
%      precMatData = zeroPrecOnePsf(PSF, center, b, tol);
%
%  Construct the data needed for a regularized circulant preconditioner
%
%  Input:
%      PSF  -  array containing a PSF
%   center  -  array containing the indices of the center of the PSF
%        b  -  blurred image to be restored
% 
%  Optional Input:
%      tol  -  truncation tolerance for preconditioner
%
%  Output:
%    precMatData  -  inverse of the (truncated) eigenvalues
%

%  J. Nagy & K. Lee  1/30/02

if nargin < 4
  tol = [];
end
center = center(1:length(size(b)));

%
% First determine if any of the dimensions of the PSF are larger than 
% the corresponding dimensions of the image, b.
% If so, we need to extract a subimage of the PSF.
%
t = fix( size(b) / 2);
t1 = max( center - t, 1 );
t2 = min( t1 + size(b) - 1, size(PSF) );
if length( size(b) ) == 1
  PSF = PSF( t1(1):t2(1) );
elseif length( size(b) ) == 2
  PSF = PSF( t1(1):t2(1), t1(2):t2(2) );
else
  PSF = PSF( t1(1):t2(1), t1(2):t2(2), t1(3):t2(3) );
end

%
% The center has moved, so let's find the new center
center = center - t1 + 1;

%
% Now determine if any of the dimensions of the PSF are smaller than
% the corresponding dimensions of the image, b.  If so, we need to
% pad the array.
%
PSF = padarray(PSF, size(b) - size(PSF), 'post');

%
% Now compute the eigenvalues of the circulant preconditioner:
%
E = circEig(PSF, center);

%
% If we were not given a truncation tolerance, let's try to
% find one:
%
if isempty(tol)
  tol = 0;  
end

if ischar(tol)
  if length(tol) > 4
    tol= defaultTol(E,b);
  else
    tol = defaultTol2(E,b);
  end
end

maxE = max( abs( E(:) ) );
E = E / maxE;                  % scale the eigenvalues, so max is 1
idx = abs(E) < tol;            % then compare to tol
E = E * maxE;                  % scale back
E(idx) = 1;                    % replace small eigenvalues by 1
E = 1 ./ E;                    % invert

precMatData = E;
