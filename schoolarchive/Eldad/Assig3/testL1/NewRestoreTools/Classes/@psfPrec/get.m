function val = get(P, prop_name)
%
%  GET Get psfPrec properties from specified object and return
%      the value.
%
%     val = get(P, prop_name);
%
%  Valid choices for prop_name are:
%    'matdata', 'type', 'boundary', 'transpose'
%

%  J. Nagy  7/8/01

switch prop_name
case 'matdata'
  val = P.matdata;
case 'type'
  val = P.type;
case 'boundary'
  val = P.boundary;
case 'transpose'
  val = P.transpose;
otherwise
  error([prop_name, 'Is not valid psfPrec property'])
end