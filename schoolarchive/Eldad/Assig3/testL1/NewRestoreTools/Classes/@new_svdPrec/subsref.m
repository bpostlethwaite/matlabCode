function B = subsref(P, index)
%
%  Define field name indexing for svdPrec object.
%
%               B = subsref(P, index)
%
%          This is called whenever a subscripted reference to the
%          svdPrec object is made, such as:
%              P(i), i = 1, 2, 3, 4
%              P.fieldname, fieldname = u, s, v
%

%  J. Nagy  7/8/01

switch index.type
case '()'
  switch index.subs{:}
  case 1
    B = P.u;
  case 2
    B = P.s;
  case 3
    B = P.v;
  otherwise
    error('Index out of range.')
  end
  
case '.'
  switch index.subs
  case 'u'
    B = P.u;
  case 's'
    B = P.s;
  case 'v'
    B = P.v;
  otherwise
    error('Invalid field names.');
  end

case '{}'
  error('Cell array indexing not supported for svdPrec object.')
end
