function [tout,yout,k1out] = ode4(ODEFUN, TSPAN, Y0, h)
%Fixed step 4th order Runge-Kutta method
%function [tout,yout] = ode4(ODEFUN, TSPAN, Y0,h)
%       ODEFUN - the function to integrate (the f in the equation dy/dx = f(x))
%       TSPAN - The time span of the integration. Can be [Tstart,Tend] or a
%               vector of timepoints to be evaluated
%       Y0    - Initial conditions (shoulve be given as y(x = 0)) 
%       h     - is the constant step

% if ~all(size(TSPAN)==[1,2])
%     
% end
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
yout(1,:) = y';

for i= 2:length(tout)
    h = tout(i)-tout(i-1);
    k1 = ODEFUN(tout(i),y);
    k2 = ODEFUN(tout(i)+h/2, y+h*k1/2);
    k3 = ODEFUN(tout(i)+h/2, y+h*k2/2);
    k4 = ODEFUN(tout(i)+h, y+h*k3);
    y = y + h*(k1 + 2*k2 + 2*k3 + k4)/6';
    yout(i,:) = y';
    k1out(i,:) = k1;
end