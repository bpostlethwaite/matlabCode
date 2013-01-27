function y = mtimes(arg1, arg2)
%
%  Overload backslash operation for svdPrec object.
%
%  Implements P * x and P' * x  for svdPrec object
%  P and vector x.  Result is returned as a vector y. 
%

%  J. Nagy  6/22/02

if ( isa(arg1, 'svdPrec') )
  [m, n] = size(arg2);
  s = reshape(arg1.s, m, n);
  y = arg1.u * ( (arg1.v' * arg2) .* s );

else

  error('incorrect argument type')

end