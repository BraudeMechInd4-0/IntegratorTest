function [tout,xout,error] = odeMPCI(ODEFUN, TSPAN, X0,AbsTol,RelTol,N,xtinit)
%[t,x,err] = odeMPCI(ODEFUN, TSPAN, X0,AbsTol,RelTol,N,xtinit)
%Bai's [1] MPCI integration method, as described in [2]
% Useage
%       ODEFUN - the function to integrate (the f in the equation dx/dt = f(x,t))
%       TSPAN  - The time span of the integration. Can be [Tstart,Tend] or a
%                vector of timepoints to be evaluated
%       X0     - Initial conditions (shoulve be given as x(t = 0))
%       AbsTol - Absolute tolerance
%       RelTol - Relative tolerace
%       N      - number of timepoints to be sampled. Optional if TSPAN is a
%                 vector longer than 2, mandatory if it is [Tstart,Tend]
%       xtinit - An initial guess. For astrodynamics it is customery to use
%                the solution for the two body problem, or even better one
%                of the SGP analytical methods.
%       tout   - A list of time points
%       xout   - Values of the function approximation on these points
%       error  - error code: 0 if all went well. -1 if failed to converge
%
%[1] Bai, X. and Junkins, J.L., 2011. Modified Chebyshev-Picard iteration methods for solution of boundary value problems. The Journal of the astronautical sciences, 58(4), pp.615-642.
%[2] Woollands, R. and Junkins, J.L., 2019. Nonlinear differential equation solvers via adaptive Picard–Chebyshev iteration: Applications in astrodynamics. Journal of Guidance, Control, and Dynamics, 42(5), pp.1007-1022.
%See also:  Darin Koblick (2024). Vectorized Picard-Chebyshev Method (https://www.mathworks.com/matlabcentral/fileexchange/36940-vectorized-picard-chebyshev-method), MATLAB Central File Exchange. Retrieved October, 2023. 

%% input check and initialization
error = 0;
if nargin <3
    error("not enough input arguments");
end
if nargin <=3
    AbsTol = 1e-12;
    RelTol = 1e-9;
end
if nargin <= 4
   RelTol = 1e-9; 
end
if nargin <= 5
    if length(TSPAN) > 2
        N = length(TSPAN);
    else
        N = 32;
    end
end
N = N-1; %because it's 0 - N in all places so if N is 32 and we don't deduct one we get 33 points. 

tau = fliplr(cos((0:N)*(pi/N))); %Gauss lobato nodes
if nargin <= 6
    [m,n] = size(X0);
    if m > 1 && n >=1 %this means it's a matrix, so we have different X's so X0 is not the boundery value but an initial guess.
        error("if you want an initial guess use xtinit")
    end
    if n == 1 && m ~=1 %it is a column instead of a vector.
        X0 = X0';
    end

    xtinit = repmat(X0,[length(tau),1]); %innitial guess
end

X0 = [X0;zeros(N,length(X0))];

xold = xtinit;%zeros(length(tau),length(X0)); %innitial guess
om2 = (TSPAN(end)-TSPAN(1))/2;
om1 = (TSPAN(end)+TSPAN(1))/2;

%% First Creat the Vectors and matrices
W = eye(length(tau));
W(1,1) = 0.5;
W(end,end) = 0.5;
T = cos((0:N-1)'*acos(tau))'; %eq A6 k = 0,1,2,...N-1
Tm1 = cos((0:N).*pi()); %acos(-1) = pi;
L = [Tm1;zeros(N,N+1)];
s_ = 1./(4:2:2*N);
S_3 = zeros(N+1);
S_3(2:end,2:end) = diag([-.5,-s_(1:end-1)],1);
S_2 = diag([1,s_],-1);
S_1 = S_2+S_3;
S = S_1(:,1:N);
S(1,:) = [1/4, zeros(1,N-1)];
A = (T'*W*T)\T'*W;


%% Get the function g using picard iterations

T = cos((0:N)'*acos(tau))';   
eAbs = inf;
eRel = inf;
i = 0;
while (eAbs  > AbsTol || eRel > RelTol) && i<2000
    F = ODEFUN(om2.*tau+om1,xold); % the VMPCM uses F = ode(input{:}).*omega2; because dx/dtau = dx/dt.dtau/dt, but APC seemed to have accounted for that in the next line
    P1 = om2*(eye(N+1) - L) * S;
    bi = X0 + P1 * A * F;
    xnew = T*bi;
    eAbs = max(abs(xnew - xold),[],'all');
    eRel = max(abs(xnew - xold)./min(xnew,xold),[],'all');
    xold = xnew;
    i = i+1;
end

if i >= 4000
    %error("did not converge")
    error = -1;
    tout = 0;
    xout=0;
    return
end


%% now evaluate at the required t's
%for each t do x(t) = sum(beta_k*T_k(t))

if length(TSPAN) ~= 2
    %evaluate at the points in TSPAN and return
    tau = -(TSPAN(end)-TSPAN(1))/(TSPAN(end)+TSPAN(1)) +2*(TSPAN)/(TSPAN(end)-TSPAN(1));
    tout = TSPAN;
else
    tout = TSPAN(1) + (tau+1)*(TSPAN(end)-TSPAN(1))/2;
end

T = cos((0:N)'.*acos(tau))'; %length(tau)xN+1
xout = zeros(length(tau),6);
for i = 1:length(tau)
    xout(i,:) = T(i,:)*bi;
end

