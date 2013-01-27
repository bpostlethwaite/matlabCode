function P = set(P, varargin)
%
% Set svdPrec properties and return the updated object.
%
%     P = set(P, ...)
% 
%  This funtion accepts an svdPrec object, P, and a variable list of
%  property name/value pairs and returns the modified object.
%  Valid properties are:
%    'u', 's', 'v'
%

%  J. Nagy  6/22/02

property_argin = varargin;

while length( property_argin ) >= 2
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);

  switch prop
  case 'u'
    P.u = val;
  case 's'
    P.s = val;
  case 'v'
    P.v = val;
  otherwise
    error('invalid property')
  end
end
