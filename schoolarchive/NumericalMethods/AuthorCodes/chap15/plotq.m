function plotq (mesh)
%
%  plot the mesh and integrand func
%

  np1 = length(mesh);
  for i=1:np1
    fx(i) = func(mesh(i));
    t(i) = 0;
  end

  plot (mesh,t,'+',mesh,fx)
 