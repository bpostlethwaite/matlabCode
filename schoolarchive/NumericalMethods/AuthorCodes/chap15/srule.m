function val = srule (a,b,fa,fab2, fb)
%
% function val = srule (a,b,fa,fab2, fb)
%
% evaluate basic Simpson rule
%
  val = (b-a)/6 * (fa + 4*fab2 + fb);
  