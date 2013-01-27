% Example 7.17 Figure 7.11 : Poisson on a square using CG, PCG/IC(tol) and multigrid.
% Try a few grids and watch iteration count and timing
clear all
close all

for l = 5:9
  N = 2^l - 1;
  Ns(l-4) = N;  % save for plotting
  tol = 1.e-6;  % should actually depend on N but never mind.

  A = delsq(numgrid('S',N+2));
  n = size(A,1);
  b = A*ones(n,1);

  % Solve using multigrid
  t0 = cputime;
  xmg = zeros(n,1); bb = norm(b);
  flevel = log2(N+1);
  for itermg = 1:100
    [xmg,rmg(itermg)] = poismg(A,b,xmg,flevel);
    if rmg(itermg)/bb < tol , break, end
  end
  time(1,l-4) = cputime - t0;
  itns(1,l-4) = itermg; 

  % CG straight
  t0 = cputime;
  [xCG,flagCG,relresCG,iterCG,rCG] = pcg(A,b,tol,2000);
  time(2,l-4) = cputime - t0;
  itns(2,l-4) = iterCG;

  % preconditioned CG with incomplete Cholesky
  t0 = cputime;
  R=cholinc(A,.01);
  [xCGpc,flagCGpc,relresCGpc,iterCGpc,rCGpc] = pcg(A,b,tol,2000,R',R);
  time(3,l-4) = cputime - t0;
  itns(3,l-4) = iterCGpc;
  
  % preconditioned CG with V-cycle
  t0 = cputime;
  x0 = zeros(n,1);
  [xcgmg,itercgmg,rcgmg] = pcgmg(A,b,x0,tol,100);
  time(4,l-4) = cputime - t0;
  itns(4,l-4) = itercgmg;
  
  figure(l-4)
  semilogy(rCG/norm(b),'r-.')
  hold on
  semilogy(rCGpc/norm(b),'m')
  semilogy(rcgmg/norm(b),'g--')
  semilogy(rmg/norm(b),'b')
  legend('CG','PCG/IC(0)','PCG/MG','MG')
  xlabel('Iterations')
  ylabel('Residuals')
  
end

% Plot relative residuals
  
  figure(10)
  plot(Ns,itns(2,:),'r-.')
  hold on
  plot(Ns,itns(3,:),'m')
  plot(Ns,itns(4,:),'g--')
  plot(Ns,itns(1,:),'b')
  
  legend('CG','PCG/IC(0)','PCG/MG','MG')
  ylabel('Iterations')
  xlabel('N')
  
  figure(11)
  plot(Ns,time(2,:),'r-.')
  hold on
  plot(Ns,time(3,:),'m')
  plot(Ns,time(4,:),'g--')
  plot(Ns,time(1,:),'b')
  legend('CG','PCG/IC(0)','PCG/MG','MG')
  ylabel('CPU secs')
  xlabel('N')