% Example 7.13 -- Figure 7.8 : convection-diffusion equation

% problem parameters
beta=.1; gamma=.1;

N=100; ee=ones(N,1);

% set up linear system matrix A
a=4; b=-1-gamma; c=-1-beta; d=-1+beta; e=-1+gamma;
t1=spdiags([c*ee,a*ee,d*ee],-1:1,N,N);
t2=spdiags([b*ee,zeros(N,1),e*ee],-1:1,N,N);
A=kron(speye(N),t1)+kron(t2,speye(N));

b=A*ones(10000,1);

% gmres(20)
[x,flag,relres,iter,resvec1]=gmres(A,b,20,1e-8,1000);

% ILU-preconditioned gmres(20)
[L,U]=luinc(A,0.01);
[x,flag,relres,iter,resvec2]=gmres(A,b,20,1e-8,1000,L,U);

semilogy(resvec1);
hold on
semilogy(resvec2,'--r')
legend('non-preconditioned GMRES','ILUT-preconditioned GMRES');
xlabel('Iteration')
ylabel('Relative residual')

