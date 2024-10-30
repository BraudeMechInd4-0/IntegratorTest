function [rt,vt] = propKep(r0,v0,t,mu,method,tol,maxIt)
%[rt,vt] = propKep(r0,v0,t,mu,method,tol,maxIt)
%This function propagates satellite position and velocity by solving
%kepler's equation

if nargin < 4
    error('not enough inputs');
elseif nargin < 5
    method = 2;
    tol = 1e-14;
    maxIt = 10000;
elseif nargin < 6
    tol = 1e-14;
    maxIt = 10000;
elseif nargin < 7
    maxIt = 10000;
elseif nargin ~= 7
    error('too many inputs?');
end

[a,e,i,O,w,f0] = kep_elements(r0,v0,mu);

% First compute E0 and the period
%E0 = acos((e+cos(f0))/(1+e*cos(f0))); %acos works between -pi and pi
%if f0 > pi
%    E0 = 2*pi - E0;
%end
E0 = wrapTo2Pi(2*pi+atan2(sin(f0)*sqrt(1-e^2),e+cos(f0)));
T = 2*pi*sqrt((a^3)/mu);
n = 2*pi/T;
M0 = E0-e*sin(E0);
t0 = M0/n; %making the perigee at some T = 0;
M = n*t+M0; %from the perigee to the point we wish to propagate


%Now we need to determine how many laps the satellite will make in the time
%it takes to propagate. So that we propagate just a piece of the orbit
if t > T
   k = fix((t+t0)/T);
   M = M - 2*pi*k;
end

[Et,c] = KepEq(e,M,method,tol,maxIt);
if c >= maxIt
    warning('did not converge on t=%.2f',t);
    rt = 0;
    vt = 0;
    return
else
    ft = wrapTo2Pi(2*pi + atan2(sqrt(1-e^2)*sin(Et),cos(Et)-e));
end


[rt,vt] = posnvelos(a,e,i,O,w,ft,mu);
end