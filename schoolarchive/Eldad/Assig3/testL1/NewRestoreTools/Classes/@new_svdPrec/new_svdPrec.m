function P = new_svdPrec(varargin)
%
%  CONSTRUCTOR FOR new_svdPrec (svd preconditioner) OBJECT
%
%  The new_svdPrec class is based on a structure with four fields:
%    U - unitary matrix
%    S - vector containing singular or spectral values
%    V - unitary matrix
%
%  Calling Syntax:
%       P = new_svdPrec             (returns object with empty fields)
%       P = new_svdPrec(svdPrecObj) (returns input object)
%       P = new_svdPrec(  A, b )
%       P = new_svdPrec(  A, b, tol )
%       P = new_svdPrec(  A, b, 'help')
%       P = new_svdPrec(  A, b, 'show help')
%
%    where 
%       * P   is a new_svdPrec object
%       * A   is a psfMatrix object or a multiPsfMatrix
%       * b   is the right hand side image for the system Ax=b that
%             is being preconditioned
%       * tol is a tolerance to "regularize" the preconditioner (e.g.,
%             the Hanke, Nagy, Plemmons approach)
%             If tol = 'help' is not specified, a default will be chosen using
%             the generalized cross validation method.
%             If tol = 'show help' some plots will be displayed so that
%                       the tolerance can be analyzed.
%

%  J. Nagy  6/22/02

%  Modifications:
%    01-09-03, J. Nagy
%              Can now use 'help' without showing plots.

switch nargin

case 0
  P.u = [];
  P.s = [];
  P.v = [];
  P = class(P, 'new_svdPrec');

case 1
  if ( isa( varargin{1}, 'new_svdPrec' ) )
    P = varargin{1};
  else
    error('Incorrect argument type')
  end

case 2
  if ( isa( varargin{1}, 'psfMatrix') )
    [U, S, V] = new_svd(varargin{1}, varargin{2});
    P.u = U;
    P.s = S;
    P.v = V;
    P = class(P, 'new_svdPrec');

  elseif (isa( varargin{1}, 'multiPsfMatrix') )
     error('new_svdPrec is not yet defined for multiPsfMatrix')

  else 
    error('Wrong input type')
  end    
        
case 3
  if ( isa( varargin{1}, 'psfMatrix') )
    [U, S, V] = new_svd(varargin{1}, varargin{2});
    if (ischar(varargin{3}))
      bhat = U'*varargin{2};
      if length(varargin{3}) > 4
        trun_tol = GCVforSVD(S, bhat(:));
      else
        trun_tol = GCVforSVD2(S, bhat(:));
      end
      trun_tol = trun_tol / max(S(:));
    else
      trun_tol = varargin{3};
    end
    S = S / max(S(:));
    S(S<trun_tol) = 1;
    P.u = U;
    P.s = S;
    P.v = V;
    P = class(P, 'new_svdPrec');

  elseif ( isa ( varargin{1}, 'multiPsfMatrix') )
     error('new_svdPrec is not yet defined for multiPsfMatrix')
    
  else 
    error('Wrong input type')
  end   
  
otherwise
  error('Incorrect number of input arguments.')
end


