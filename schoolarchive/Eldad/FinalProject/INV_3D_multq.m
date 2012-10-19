% Mike McMillan
% #20741039
% MATH 607d
% Project


% This code performs 3D Inversion of DC Resistivity data
% This code uses multiple source locations

% To run this code:
% 1) Choose number of cells in x,y,z direction
% 2) Choose a true conductivity model to calculate observed data
% 3) Choose a best guess reference model
% 4) Pick location to inject current source
% 5) Add a percentage of Gaussian noise to the observed data
% 6) Choose a tolerance and maximum # of iterations for Gauss-Newton loop
% 7) Choose a starting regularization parameter alpha.

clear all
close all

%kk = 1;
%fprintf('   h     ||err r1|| ||err r2|| ||un-ua|| \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forward Modelling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% counter
kk=1;

% Define dimensions of problem
%for N = [4,8,16];
N=10;
n1=N;n2=N;n3=N;
n = n1*n2*n3;
h1 = 1/n1;
h2 = 1/n2;
h3 = 1/n3;

% Define unit vector
e = ones(n,1);

% Define model for forward modelling: block in a halfspace
con = 5; % Conductivity of block
m = ones(n,1); % background conductivity 
m = reshape(m,n1,n2,n3);

% Block dimensions: (2 cells x 2 cells x 2 cells)
% Location: centre two cells in y and z. Surface two cells in x.
m(2:3,ceil(n2/2):ceil(n2/2)+1,ceil(n3/2):ceil(n3/2)+1) = con*ones(2,2,2);

%% View model (if necessary)

%{
for i = 1:4
    uplot = squeeze(m(i,:,:));
    figure(1)
    subplot(2,2,i)
    imagesc(uplot)
    caxis([0 5]);
    colorbar
    xlabel('x');ylabel('y');
    title(['model (S/m) at depth slice = ',num2str(i)]);
end
%}

%%  Define q and Q
m = m(:);

% Define numerous transmitter electrode locations q
num_lines = 3; % y direction on slices
num_dipoles = 2; % x direction on slices
n_space = 3;
numq = num_lines*num_dipoles*n_space;
qtot = zeros(n,numq);
kkk = 1;
for i = 1:num_lines;
    for j = 1:num_dipoles
        for k = 1:n_space
            q = zeros(n1,n2,n3);
            q(1,i+4,j+2) = -1/(h1*h2*h3);
            q(1,i+4,j+k+3) = 1/(h1*h2*h3);
            q = q(:);
            qtot(:,kkk) = q;
            kkk=kkk+1;
        end
    end
end

% Define matrix Q which outlines where data will be collected
% If data is collected at every surface cell:
num_surf = n2*n3;
num_data = numq*num_surf;
Q = speye(num_surf);
k = zeros(1,n1);k(1)=1;
Q = kron(Q,k);
Q = kron(speye(numq),Q); % For multiple q purposes

%% View source locations q (change column of qtot to see diff q locations)

%
for i = 1:numq;
subplot(ceil(sqrt(numq)),floor(sqrt(numq)),i)
qtotplot = reshape(qtot(:,i),n1,n2,n3);
qplot = (squeeze(qtotplot(1,:,:)));
imagesc(qplot)
end
%}

%% Calculate forward modelled potentials and Sensitivity Matrix J
% potentials u calculated at every cell centre. 

% Calculate u and J using block in a halfspace model:
[u,J] = dc3dfwd_multq(n1,n2,n3,m,qtot(:),Q,numq);

% Create data observed by adding Gaussian noise:
d_obs = Q*u;
noise_level = 0.05;
noise = noise_level*randn(length(d_obs),1).*d_obs;
d_obs_noise = d_obs + noise;

% Calculate u and J using only halfspace model: used to see if signal
% exists by plotting the differences in u between block in halfspace model
% and just halfspace model. Use plotting tools in "plot u values"
mh = ones(n,1);
[uh,Jh] = dc3dfwd_multq(n1,n2,n3,mh,qtot(:),Q,numq);

%% Define grid points (if necessary)

%{
% cell centered grid
[xc yc zc] = ndgrid(h1/2:h1:1-h1/2,h2/2:h2:1-h2/2,h3/2:h3:1-h3/2);

% nodal grid
[xn yn zn] = ndgrid(0:h1:1,0:h2:1,0:h3:1);

% staggered grid for Jx
[xdx ydx zdx] = ndgrid(0:h1:1,h2/2:h2:1-h2/2,h3/2:h3:1-h3/2);
[xdy ydy zdy] = ndgrid(h1/2:h1:1-h1/2,0:h2:1,h3/2:h3:1-h3/2);
[xdz ydz zdz] = ndgrid(h1/2:h1:1-h1/2,h2/2:h2:1-h2/2,0:h3:1);

% boundary points
xb = [xn(1,1,1),xn(end,end,end)]';
yb = [yn(1,1,1),yn(end,end,end)]';
zb = [zn(1,1,1),zn(end,end,end)]';

%}

%% Plot u values (if necessary)

%{
% first with block in halfspace model
figure(1);
subplot(121)
u = reshape(u,n1,n2,n3);
for i = 1:1;
    uplot = squeeze(u(i,:,:));
    imagesc(uplot)
    colorbar
    pause(2)
end
u = u(:);

% second with halfspace model
m = ones(n,1);
%[uh,J] = dc3dfwd(n1,n2,n3,m,q,Q);
subplot(122)
uh = reshape(uh,n1,n2,n3);
for i = 1:1;
    uplot = squeeze(uh(i,:,:));
    imagesc(uplot)
    colorbar
    pause(2)
end
uh = uh(:);
%

% Plot differences in two models
diff_u = u(:) - uh(:);
diff_u = reshape(diff_u(1001:2000),n1,n2,n3);
for i = 1:1;
    uplot = squeeze(diff_u(i,:,:));
    figure(1)
    imagesc(uplot)
    colorbar
end
%}

%% Adjoint test for J - adj_testJ should equal 0

%{
v = randn(n*numq,1);
z = randn(num_surf*numq,1);
t1 = z'*J*v;
t2 = v'*J'*z;
adj_testJ = norm(t1-t2);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regularization - use L2 norm approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct Regularization Gradient operators in x,y,z

d1 = [ones((n1-1),1);0];
d2 = [ones((n2-1),1);0];
d3 = [ones((n3-1),1);0];
RDc2n_1 = 1/h1 * spdiags([-d1,[0;d1(1:end-1)]],(-1:0),n1+1,n1);
RDc2n_2 = 1/h2 * spdiags([-d2,[0;d2(1:end-1)]],(-1:0),n2+1,n2);
RDc2n_3 = 1/h3 * spdiags([-d3,[0;d3(1:end-1)]],(-1:0),n3+1,n3);

% Extend x,y,z gradient operators to 3D using kroenicker products
RG1 = kron(speye(n3),kron(speye(n2),RDc2n_1));
RG2 = kron(speye(n3),kron(RDc2n_2,speye(n1)));
RG3 = kron(RDc2n_3,kron(speye(n2),speye(n1)));

% Put x,y,z operators together to form Regularization Gradient operator
G = [RG1;RG2;RG3];

%% Construct Regularization Averaging Operator

RAn2c_1 = spdiags([ones(n1,1),ones(n1,1)],0:1,n1,n1+1);
RAn2c_2 = spdiags([ones(n2,1),ones(n2,1)],0:1,n2,n2+1);
RAn2c_3 = spdiags([ones(n3,1),ones(n3,1)],0:1,n3,n3+1);

% Extend to 3D using Kroenicker products
RA1 = kron(speye(n3),kron(speye(n2),RAn2c_1));
RA2 = kron(speye(n3),kron(RAn2c_2,speye(n1)));
RA3 = kron(RAn2c_3,kron(speye(n2),speye(n1)));
RAv = [RA1,RA2,RA3];

% Regularization Operator
Rm = @(model)(0.5*(h1*h2*h3)*e'*RAv*((G*model).^2));

%% Test Regularization Discretization

%{
% For this test, must turn on for loop at beginning for N = [4,8,16];
% Place an "end" at the end of the code

reg_m = @(x,y,z)(cos(2*pi.*x).*cos(2*pi.*y).*cos(2*pi.*z));
R_test = R(reg_m(xc(:),yc(:),zc(:)));
R_diff(kk) = 3*pi^2/2 - R_test;
kk = kk+1;
h=[1/4,1/8,1/16];
count = (1:3);
semilogy(count,R_diff,count,h.^2)
xlabel('Count');ylabel('Regularization Difference');
title('Order h^2 convergence of regularization operator');
legend('Regularization Differnce','Order h^2');
%}
    
%% Create grad_Rm and grad2_Rm operators
Av_e = RAv'*(h1*h2*h3)*e;
[sx,sy] = size(Av_e);
grad_Rm = @(model)(G'*(spdiags((Av_e),0,sx,sx))*(G*model));
grad2_Rm = G'*(spdiags((Av_e),0,sx,sx))*G;

%% Gauss-Newton Step:

% Perform Gauss-Newton Loop
s = 1;
tol = 1e-4;
iter = 1;
alpha = 0.01;
m_guess = mh;
m_guess = m_guess(:);

% While the model step s is above a tolerance continue iterating until the
% maximum number of iterations has been reached.

% maxiter = 6;
% stot = zeros(maxiter,1);
% while norm(s,'inf') > tol && iter < maxiter;
%     [u,J] = dc3dfwd_multq(n1,n2,n3,m_guess,qtot(:),Q,numq);
%     lhs = J'*J + alpha*grad2_Rm;
%     rhs = -J'*(Q*u - d_obs_noise) + alpha*grad_Rm(m_guess-mh);
%     s = lhs \ rhs;
%     stot(iter) = norm(s,'inf');
%     m_guess = m_guess + s;
%     alpha = alpha / 2;
%     iter = iter + 1;
% end

% Test to find Tikhonov curve - allocate memory
%alpha = logspace(2,-4,20);
misfit_m = zeros(length(alpha),1);
misfit_d = zeros(length(alpha),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test stopping criteria
%for j = 1:length(alpha)  
    m_guess = mh;
    s = 1;
    iter = 1;
    maxiter = 6;
    %mis_d = 1e9;
    fprintf('data_misfit   target_misfit  model_misfit  alpha \n')
    while iter < maxiter;% && mis_d > num_data;%norm(s,'inf') > tol && iter < maxiter;
        [u,J] = dc3dfwd_multq(n1,n2,n3,m_guess,qtot(:),Q,numq);
        misfit_m = (Rm(m_guess-mh))^2;
        misfit_d_temp = (Q*u - d_obs_noise);
        misfit_d = norm((misfit_d_temp./(noise_level*d_obs)).^2);
        %mis_d = norm((misfit_d_temp./(noise_level*d_obs)).^2);
        lhs = J'*J + alpha*grad2_Rm;
        rhs = -J'*misfit_d_temp + alpha*grad_Rm(m_guess-mh);
        s = lhs \ rhs;
        m_guess = m_guess + s;
        iter = iter + 1;
        norm(s,'inf');
        fprintf('%e %e %e %e \n',misfit_d,length(d_obs),misfit_m,alpha)
        alpha = alpha / 2;
    end
%end

figure(9)
plot(misfit_m,misfit_d,'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot results
m = reshape(m,n1,n2,n3);
m_guess = reshape(m_guess,n1,n2,n3);
for i = 1:4
    uplot = squeeze(m_guess(i,:,:));
    figure(3)
    subplot(2,2,i)
    imagesc(uplot)
    caxis([0 5]);
    colorbar
    xlabel('x');ylabel('y');
    title(['model (S/m) at depth slice = ',num2str(i)]);
end


