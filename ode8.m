function [tout,yout,k1out] = ode8(ODEFUN, TSPAN, Y0, h)
%Fixed step 8th order Runge-Kutta method
%function [tout,yout] = ode8(ODEFUN, TSPAN, Y0,h)
%       ODEFUN - the function to integrate (the f in the equation dy/dx = f(x))
%       TSPAN - The time span of the integration. Can be [Tstart,Tend] or a
%               vector of timepoints to be evaluated
%       Y0    - Initial conditions (shoulve be given as y(x = 0)) 
%       h     - is the constant step

c = [0 1/9 1/6  1/4 1/10 1/6 1/2 2/3 1/3 5/6 5/6 1];
b = [41/840 0 0 0 0 216/840 272/840 27/840 27/840 36/840 180/840 41/840];
a = [zeros(1,12);...
    1/9, zeros(1,11);...
    1/24, 3/24, zeros(1,10);...
    1/16, 0, 3/16, zeros(1,9);...
    29/500, 0, 33/500, -12/500,zeros(1,8);...
    33/972, 0, 0, 4/972, 125/972, zeros(1,7);...
    -21/36, 0, 0, 76/36, 125/36, -162/36, zeros(1,6);...
    -30/243, 0, 0, -32/243, 125/243, 0, 99/243, zeros(1,5);...
    1175/324, 0, 0, -3456/324, -6250/324, 8424/324, 242/324, -27/324, zeros(1,4);...
    293/324, 0, 0, -852/324, -1375/324, 1836/324, -118/324, 162/324, 1, zeros(1,3);...
    1303/1620, 0, 0, -4260/1620, -6875/1620, 9990/1620, 1030/1620, 0, 0, 162/1620, 0, 0;...
    -8595/4428, 0, 0, 30720/4428, 48750/4428, -66096/4428, 378/4428, -729/4428, -1944/4428, -1296/4428, 3240/4428, 0];


if nargin > 3 && length(TSPAN) == 2
    tout = TSPAN(1):h:TSPAN(2)';
elseif length(TSPAN) < 3
    error("size of TSPAN wrong, should be [T0 TFINAL] or contain many points");
else
    tout = TSPAN;
    if size(tout,1) > 1
        tout = tout';
    end   
end
yout = zeros(length(tout),length(Y0));
k1out = zeros(length(tout),length(Y0));

y = Y0';
yout(1,:) = Y0; %Y0 is 1x6
%k1out(1,:) = ODEFUN(tout(1),Y0); %Y0 is 1x6

for i= 2:length(tout)
    h = tout(i)-tout(i-1);
    %computke the k's
    k = ODEFUN(tout(i),y); %k1 6x1
    for s = 2:length(b)
        tmp = sum(a(s,1:s-1).*k,2);
        k(:,s) = ODEFUN(tout(i)+c(s)*h,y+h*tmp);
    end
    %sum up for the next yn+1
    phi = sum(b.*k,2);
    y = y + h*phi;
    yout(i,:) = y';
    k1out(i-1,:) = k(:,1)';
end