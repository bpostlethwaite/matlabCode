function A=multiPsfMatrix(varargin)
%
%  CONSTRUCTOR FOR multiPsfMatrix OBJECT
%
%  The multiPsfMatrix class is based on a structure with three fields:
%
%    psfMatrices - cell array containing a bunch of psfMatrix objects
%    transpose   - indicates if the matrix has been transposed.
%    boundary    - character string array indicating type of boundary conditions
%                  to be used.  Choices are:
%                  'zero', 'periodic', 'reflexive' (or 'neumann')
%                  The default is zero.
%
%  Calling Syntax:
%       A = multiPsfMatrix             (returns object with empty fields)
%       A = multiPsfMatrix(multiPsfMatrixObj) 
%       A = psfMatrix(K1,K2, ...)
%       A = psfMatrix(A,K1,K2, ...)
%
%    where 
%       * multiPsfMatrixObj is an already existing multiPsfMatrix object
%       * K1, K2, ...       are all psfMatrix objects
%       * A                 is a multiPsfMatrix object
%

%  J. Nagy & K. Lee  2/8/02

if nargin == 1
  if (isa(varargin{1}, 'cell') )
    varargin = varargin{1};
    nn = length(varargin);
  else 
    nn = 1;
  end
else
  nn = length(varargin);
end

switch nn

case 0
  A.psfMatrices = cell(0);
  A.transpose = 0;
  %A.boundary = '';
  A = class(A, 'multiPsfMatrix');

case 1
  if ( isa(varargin{1}, 'multiPsfMatrix') )
    A=varargin{1};
  elseif ( isa( varargin{1}, 'psfMatrix') )
    c = cell(1);
    c{1} =varargin{1};
   % d = varargin{1};
    A.psfMatrices = c;
    A.transpose = 0;
    %A.boundary = d.boundary;
   % A= class(A, 'multiPsfMatrix');
  elseif ( isa( varargin{1}, 'double') )
    c{1} = psfMatrix( varargin{1} );
   % d=c{1};
    A.psfMatrices = c;
    A.transpose = 0;
    %d=c{1};
    %A.boundary=d.boundary;
   % A = class(A, 'multiPsfMatrix');
  else
    error('**Incorrect argument type')
  end
     A = class(A, 'multiPsfMatrix');

otherwise,
  n = length(varargin);
  if ( isa( varargin{1}, 'multiPsfMatrix') )
    A1 = varargin{1};
    %b=A1.boundary
    L1 = length(A1.psfMatrices);
    c = cell(L1+n-1,1);
    c(1:L1) = A1.psfMatrices;
    for i = 2:n
      if ( isa( varargin{i},'psfMatrix') )
        c{L1+i-1} = varargin{i};
        %d=c{L1+i-1};
        %b1=d.boundary;
        %if b ~= b1
        %    error('boundaries must be the same');
        %end
      else
        error('****Incorrect argument type')
      end
    end
  elseif ( isa( varargin{1},'double' ))
    c = cell(n,1);
    for i = 1:n
      if (isa (varargin{i},'double'))
        c{i} = psfMatrix( varargin{i} );
        %b=c{1};
        %d=c{i};
       % b1=d.boundary;
       % if b ~= b1
         %   error('boundaries must be the same');
         %end
      else
        error('******Incorrect argument type')
      end
    end
  else
    c = cell(n,1);
   % d=varargin{1};
   % b=d.boundary;
    for i = 1:n
      if (isa (varargin{i},'psfMatrix'))
        c{i} = varargin{i};
        %d=c{i};
       % b1=d.boundary;
       % if b ~= b1
        %    error('boundaries must be the same');
        %end
      else
        error('********Incorrect argument type')
      end
    end
  end
  A.psfMatrices = c;
  A.transpose = 0;
 % A.boundary = b;
  A = class(A, 'multiPsfMatrix');
end
