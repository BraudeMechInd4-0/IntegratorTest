function [En, counter] = KepEq(e,M,method,tol,maxIt)
%[E] = KepEq(e,M) or [En, counter] = KepEq(e,M,tol,maxIt)
%Solve keplers Eq.
%method - initial start point, should be 1: normal, 2: fast, 3: faster? 31, 21
%variations that are faster in some cases
%e - eccentricity
%M - mean anomaly
%E - eccentric anomaly
%E = M+esin(E)

if nargin < 4
    tol = 1e-14;
    maxIt = 10000;
elseif nargin < 3
    tol = 1e-14;
    maxIt = 10000;
    method = 2;
elseif nargin ~=5
    error("wrong number of inputs");
end

%counter = ones(1,length(M));

err = 1000;
switch method
    case 1
        En = M;
    case 2
        En = M+e*sin(M)+0.5*(e^2)*sin(2*M);
    case 3
        En = M + e*sin(M) + 0.5*(e^2)*sin(2*M)+(3/8*sin(3*M)-1/8*sin(M))*e^3;
    case 21
        En = M+e*cos(M)+0.5*(e^2)*cos(2*M);
    case 31
        En = M+e*cos(M)+0.5*(e^2)*cos(2*M) + 1/(2*3)*(e^3)*cos(3*M);
    otherwise
        error("bad method");
end


counter = 1;


while err > tol && counter < maxIt
    Enp1 = En + (M - En + e*sin(En))/(1 - e*cos(En));
    err = abs(Enp1 - En);
    En = Enp1;
    counter = counter + 1;
end
if ~isnumeric(err)
    counter = maxIt+1;
end
