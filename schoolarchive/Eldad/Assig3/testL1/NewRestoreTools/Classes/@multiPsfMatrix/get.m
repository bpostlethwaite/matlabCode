function val = get(A, prop_name)
%
%  Get psfMatrix properties from specified object and return
%  the value.
%
%     val = get(A, prop_name);
%
%  Valid choices for prop_name are:
%    'psfMatrices', 'boundary', 'transpose'
%

%  K. Lee 2/6/02

switch prop_name
case 'psfMatrices'
  val = A.psfMatrices;
case 'boundary'
  val = A.boundary;
case 'transpose'
  val = A.transpose;
otherwise
  error([prop_name, 'Is not valid psfMatrix property'])
end