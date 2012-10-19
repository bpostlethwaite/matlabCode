function val = get(P, prop_name)
%
%  GET Get svdPrec properties from specified object and return
%      the value.
%
%     val = get(P, prop_name);
%
%  Valid choices for prop_name are:
%    'u', 's', 'v'
%

%  J. Nagy  6/22/02

switch prop_name
case 'u'
  val = P.u;
case 's'
  val = P.s;
case 'v'
  val = P.v;
otherwise
  error([prop_name, 'Is not valid svdPrec property'])
end