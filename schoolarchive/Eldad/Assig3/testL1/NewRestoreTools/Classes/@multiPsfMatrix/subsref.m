function B = subsref(A, index)
%
%  Define field name indexing for psfMatrix object.
%
%               B = subsref(A, index)
%
%          This is called whenever a subscripted reference to the
%          psfMatrix object is made, such as:
%              A(i)
%              A.fieldname, fieldname = psfMatrices, boundary, transpose
%

%  K. Lee 2/8/02

switch index.type
case '()'
  switch index.subs{:}
  case 1
    B = A.psfMatrices;
    % case 2
  %  B = A.boundary;
  case 2
    B = A.transpose;
  otherwise
    error('Index out of range.')
  end
  
case '.'
  switch index.subs
  case 'psfMatrices'
    B = A.psfMatrices;
    %case 'boundary'
  %  B = A.boundary;
  case 'transpose'
    B = A.transpose;
  otherwise
      error('blah')
   % B = subsref(A.psfMatrices, index);
  end

case '{}'
  error('Cell array indexing not supported for psfMatrix object.')
end
