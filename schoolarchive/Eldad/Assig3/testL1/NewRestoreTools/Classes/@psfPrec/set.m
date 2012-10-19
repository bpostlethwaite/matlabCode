function P = set(P, varargin)
%
% Set psfPrec properties and return the updated object.
%
%     P = set(P, ...)
% 
%  This funtion accepts an psfPrec object, P, and a variable list of
%  property name/value pairs and returns the modified object.
%  Valid properties are:
%    'matdata', 'type', 'boundary', 'transpose'
%

%  J. Nagy  7/8/01

property_argin = varargin;

while length( property_argin ) >= 2
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);

  switch prop
  case 'matdata'
    P.matdata = val;
  case 'type'
    P.type = val;
  case 'boundary'
    P.boundary = val;
  case 'transpose'
    P.transpse = val;
  otherwise
    error('invalid property')
  end
end
