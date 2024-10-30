function u = hermiteInterp(y0,y1,f0,f1,theta,h)
% function u = hermiteInterp(y0,y1,f0,f1,theta,h)
%For an RK method, iterpolate using known points
% Usage:
%       y0,y1 -  given approximation on two points
%       f0,f1 -  the first derivatives on these points
%       h     -  the step between the time points.
%       thata -  the point to which we wish to interpolate, t =
%       t_0+theta*h
%

u = (1-theta)*y0+theta*y1+theta*(theta-1)*((1-2*theta)*(y1-y0)+(theta-1)*h*f0+theta*h*f1);
