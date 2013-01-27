% Example 15.2 simple basic rule integration of f(x) = exp(x) 

case1 = 2;

if case1 == 1
  % over [0,1]
  f1 = exp(1); fh = exp(.5); f0 = 1;

  Iexact = f1 - f0;
  Itrap = .5*(f0 + f1)
  Etrap = Iexact - Itrap
  Imid = fh
  Emid = Iexact - Imid
  Isimp = (f0 + 4*fh + f1) / 6
  Esimp = Iexact - Isimp

else

  % over [0.9,1]
  f1 = exp(1); fh = exp(.95); f0 = exp(.9);

  Iexact = f1 - f0;
  Itrap = .5*(f0 + f1) * .1
  Etrap = Iexact - Itrap
  Imid = fh * .1
  Emid = Iexact - Imid
  Isimp = (f0 + 4*fh + f1) / 6 * .1
  Esimp = Iexact - Isimp
  
end