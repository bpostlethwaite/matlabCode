function y = mldivide(arg1, arg2)
%
%  Overload backslash operation for svdPrec object.
%
%  Implements P \ b and P' \ b   for svdPrec object
%  P and vector b.  Result is returned as a vector y. 
%

%  J. Nagy  6/22/02

if ( isa(arg1, 'svdPrec') )
  %  This is a stupid hack to make this function work for
  %  vector inputs for arg2.
  [m, n] = size(arg2);

  if n == 1
    U = arg1.u;
    Ua = U.a;
    Ub = U.b;
    mm = size(Ub, 1);
    nn = size(Ua, 1);
    arg2 = reshape(arg2, mm, nn);
    m = mm;, n = nn;
    rs = true;
  else
    rs = false;
  end
  s = reshape(arg1.s, m, n);
  y = arg1.v * ( (arg1.u' * arg2) ./ s );

else

  error('incorrect argument type')

end
if rs
  y = y(:);
end
