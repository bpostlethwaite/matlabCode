function B = subsref(P, index)
%
%  Define field name indexing for psfPrec object.
%
%               B = subsref(P, index)
%
%          This is called whenever a subscripted reference to the
%          psfPrec object is made, such as:
%              P(i), i = 1, 2, 3, 4
%              P.fieldname, fieldname = matdata, type, boundary, transpose
%

%  J. Nagy  7/8/01

switch index.type
case '()'
  switch index.subs{:}
  case 1
    B = P.matdata;
  case 2
    B = P.type;
  case 3
    B = P.boundary;
  case 4
    B = P.transpose;
  otherwise
    error('Index out of range.')
  end
  
case '.'
  switch index.subs
  case 'matdata'
    B = P.matdata;
  case 'type'
    B = P.type;
  case 'boundary'
    B = P.boundary;
  case 'transpose'
    B = P.transpose;
  otherwise
    error('Invalid field names.');
  end

case '{}'
  error('Cell array indexing not supported for psfPrec object.')
end
