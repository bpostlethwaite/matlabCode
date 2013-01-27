function [x_out, exit] = HyBR(A, b, P, options)
%
% [x_out, exit] = HyBR(A, b, P, options)
%
% HyBR is a Hybrid Bidiagonalization Regularization method used for 
% solving large-scale, ill-posed inverse problems of the form:
%               b = A*x + noise
% The method combines an iterative Lanczos Bidiagonalization (LBD) Method 
% with a SVD-based regularization method to stabilize the semiconvergence
% behavior that is characteristic of many ill-posed problems.
%
% Inputs:
%       A : either (a) a full or sparse matrix
%                  (b) a matrix object that performs matrix*vector and
%                       matrix'*vector operations
%       b : rhs vector or vector object
%       P : preconditioner (optional)
% options : structure with the following fields (optional)
%         InSolv - solver for the inner problem: [none | {Tikhonov}]
%         RegPar - specify a value, or a method to find reg. parameter:
%                       [value | GCV | {modGCV} | optimal]
%                   Note: 'optimal' requires x_true
%          Omega - if RegPar is 'modGCV', then omega must be 
%                       [value | {adapt}]
%           Iter - maximum number of Lanczos iterations: 
%                       [ value | {min(m,n,100)} ]
%         Reorth - reorthogonalize Lanczos subspaces: [on | {off}]
%         x_true - True solution : [ array | {off} ]
%                Returns error with respect to x_true at each iteration
%                and is used to compute 'optimal' regularization parameters
%         BegReg - Begin regularization after this iteration: 
%                   [ positive scalar | {2} ]
%             Vx - extra space needed for finding optimal reg. parameters
%     
%       Note: options is a structure created using the function 'HyBRset' 
%               (see 'HyBRset' for more details)
%
% Outputs:
%      x_out : computed solution
%       exit : structure with the following fields:
%      iterations - stopping iteration (options.Iter | GCV-determined)
%         GCVstop - GCV curve used to find stopping iteration
%            Enrm - relative error norms (requires x_true)
%            Rnrm - relative residual norms
%            Xnrm - relative solution norms
%             U,V - Lanczos basis vectors
%               B - bidiagonal matrix from LBD
%            flag - a flag that describes the exit/stopping condition:
%                       1 - flat GCV curve 
%                       2 - min of GCV curve (within window of 4 its)
%                       3 - performed max number of iterations
%
% References:
%   [1] Paige and Saunders, "LSQR an algorithm for sparse linear
%       equations an sparse least squares", ACM Trans. Math Software,
%       8 (1982), pp. 43-71.
%   [2] Bjorck, Grimme and Van Dooren, "An implicit shift bidiagonalization
%       algorithm for ill-posed systems", BIT 34 (11994), pp. 520-534.
%   [3] Chung, Nagy and O'Leary, "A Modified GCV Method for Lanczos Hybrid
%       Regularization", in process.
%

%% Initialization
defaultopt = struct('InSolv','tikhonov','RegPar','modgcv','Omega',...
  'adapt', 'Iter', [] , 'Reorth', 'off', 'x_true', 'off', 'BegReg', 2,...
  'Vx' , []);

% If just 'defaults' is passed in, return the default options in x_out
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    x_out = defaultopt;
    return;
end

% Check for acceptable number of input arguments
if nargin < 2
  error('HyBR: Not Enough Inputs')
elseif nargin < 3
  P = []; options = [];
elseif nargin < 4
  options = [];
end
if isempty(options)
  options = defaultopt;
end

% Get options:
[m,n] = size(A); defaultopt.Iter = min([m, n, 100]);
options = HyBRset(defaultopt, options);
solver = HyBRget(options,'InSolv',[],'fast');
regpar = HyBRget(options,'RegPar',[],'fast');
omega = HyBRget(options,'Omega',[],'fast');
maxiter = HyBRget(options,'Iter',[],'fast');
x_true = HyBRget(options,'x_true',[],'fast');
regstart = HyBRget(options,'BegReg',[],'fast');

adaptmodGCV = strcmp(regpar, {'modgcv'}) && strcmp(omega, {'adapt'});
notrue = strcmp(x_true,{'off'});
if ~notrue
  nrmtrue = norm(x_true);
end

% Set-up exit parameters:
exitparams = nargout>1;
if exitparams 
  exit.iterations = maxiter;
  exit.GCVstop = [];
  exit.Enrm = ones(maxiter,1);
  exit.Rnrm = ones(maxiter,1);
  exit.Xnrm = ones(maxiter,1);
  exit.U = [];
  exit.V = [];
  exit.B = [];
  exit.flag = 3;
end

% Test for a preconditioner:
if isempty(P)
  beta = norm(b); U = b / beta;
  handle = @LBD;
else
  U = P\b;
  beta = norm(U); U = U / beta;
  handle = @PLBD;
end

%% Main Code
B = []; V = []; GCV = []; Omega= []; x_out = [];
insolve = 'none'; terminate = 1; warning = 0;

h = waitbar(0, 'Beginning iterations: please wait ...');

for i = 1:maxiter+1 %First iteration (i=1) is just initialization
  [U, B, V] = feval(handle, A, U, B, V, P, options);
  vector = (beta(1)*eye(size(U,2),1));
  
  if ~notrue
    options.Vx = V'*x_true;
  end

  if i >= 2 %Begin Lanczos iterations
    if i >= regstart %Begin to regularize projected problem now
      insolve = solver;
    end
    switch insolve
      case 'tikhonov'
        [Ub, Sb, Vb] = svd(B);

        if adaptmodGCV %if we use the adaptive, modified GCV method
          Omega(i-1) = min(1, findomega(Ub'*vector, diag(Sb)));
          options.Omega = mean(Omega);
        end
        
        % Solve the projected problem with Tikhonov
        [f, alpha] = Tikhonovsolver(Ub, diag(Sb), Vb, vector, options);

        % Compute the GCV value used to find the stopping criteria
        [GCV(i-1),Uend(i-1)] = GCVstopfun(alpha, Ub(1,:)', diag(Sb), beta, n);
        
        % Determine if GCV wants us to stop
        if i>2 && terminate
          %%-------- GCV curve is flat, we stop ------------------------
          if abs((GCV(i-1)-GCV(i-2)))/GCV(1) < 10^-4  
            x_out = V*f; 
            if notrue %Set all the exit parameters and return
              if exitparams 
                exit.U = U;
                exit.V = V;
                exit.B = B;
                exit.GCVstop = GCV(:);
                exit.iterations = i-1;
                exit.flag = 1;
              end
              close(h)
              return;
            else % Flat GCV curve means stop, but continue since have x_true
                if exitparams
                  exit.iterations = i-1;
                  exit.flag = 1;
                end            
            end
            terminate = 0;

          %%--- Have warning: Avoid bumps by using a window of 4 its --
          elseif warning && length(GCV) > iterations_save + 3
            if GCV(iterations_save) < GCV(iterations_save+1:end)
              % We should have stopped at iterations_save.
              x_out = x_save; 
              if notrue %Set all the exit parameters and return
                if exitparams 
                  exit.U = U;
                  exit.V = V;
                  exit.B = B;
                  exit.GCVstop = GCV(:);
                  exit.iterations = iterations_save;
                  exit.flag = 2;
                end
                close(h)
                return;
              else % GCV says should stop at iterations_save, but continue since have x_true
                if exitparams
                  exit.iterations = iterations_save;
                  exit.flag = 2;
                end
              end
              terminate = 0; %Solution is already found!

            else % It was just a bump... keep going
              warning = 0;
              x_out = [];
              iterations_save = maxiter;
            end
          %% ----- No warning: Check GCV function-----------------------
          elseif ~warning 
            if GCV(i-2) < GCV(i-1)
              warning = 1;
              x_save = V*f;
              iterations_save = i-1;
            end
          end
        end
 
      case 'none'
        f = B \ vector;    
        
      otherwise
        error('HyBR error: No inner solver!')
    end
    x = V*f; 
    
    if exitparams
      if ~notrue
        exit.Enrm(i-1,1) = norm(x-x_true)/nrmtrue;
      end
      exit.Rnrm(i-1,1) = norm(b(:) - A*x(:));
      exit.Xnrm(i-1,1) = norm(x(:));
    end
    
  end
  waitbar(i/(maxiter+1), h)
end
close(h)

if isempty(x_out) %GCV never terminated the process
  x_out = x;
end

if exitparams
  exit.U = U;
  exit.V = V;
  exit.B = B;
  exit.GCVstop = GCV(:);
end



%% -----------------------SUBFUNCTION---------------------------------------
function omega = findomega(bhat, s)
%
%   omega = findomega(bhat, s)
%
%  This function computes a value for the omega parameter.
%
%  The method: Assume the 'optimal' regularization parameter to be the
%  smallest singular value.  Then we take the derivative of the GCV
%  function with respect to alpha, evaluate it at alpha_opt, set the 
%  derivative equal to zero and then solve for omega.
%  
%  Input:   bhat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%
%  Output:     omega - computed value for the omega parameter.

%
%   First assume the 'optimal' regularization parameter to be the smallest
%   singular value.
%
alpha = s(end);

%
% Compute the needed elements for the function.
%
m = length(bhat);
n = length(s);
t0 = sum(abs(bhat(n+1:m)).^2);

s2 = abs(s) .^ 2;
alpha2 = alpha^2;

tt = 1 ./ (s2 + alpha2);

t1 = sum(s2 .* tt);
t2 = abs(bhat(1:n).*alpha.*s) .^2;
t3 = sum(t2 .* abs((tt.^3)));

t4 = sum((s.*tt) .^2);
t5 = sum((abs(alpha2*bhat(1:n).*tt)).^2);

v1 = abs(bhat(1:n).*s).^2;
v2 = sum(v1.* abs((tt.^3)));

%
% Now compute the omega value.
%
omega = (m*alpha2*v2)/(t1*t3 + t4*(t5 + t0));


%% ---------------SUBFUNCTION ---------------------------------------
function [G,uend] = GCVstopfun(alpha, u, s, beta, n)
%   
%  Evaluate the GCV function G(i, alpha) to determine a stopping index.
%
uend = u(end);
k = length(s);
beta2 = beta^2;

s2 = abs(s) .^ 2;
alpha2 = alpha^2;

t1 = 1 ./ (s2 + alpha2);
t2 = abs(alpha2*u(1:k) .* t1) .^2;
t3 = s2 .* t1;

num = beta2*(sum(t2) + abs(u(k+1))^2)/n;
den = ( (n - sum(t3))/n )^2;
G = num / den;