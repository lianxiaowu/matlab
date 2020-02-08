function[nx,tau] = quad1d(x,ss,delta)
%QUAD1D	1D quad. zero finder
%
% [nx,tau] = quad1d(x,ss,delta) tau is min(1,step to zero
% of a 1-D quadratic ay^2 + b*y + c , where 
% a = x'*x; b = 2*(ss'*x); c = ss'*ss-delta^2). nx is the
% new x value, nx = tau*x;

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1998/03/28 19:59:40 $

a = x'*x;
b = 2*(ss'*x); 
c = ss'*ss-delta^2;
numer = -(b + sign(b)*sqrt(b^2-4*a*c));
if a == 0 & (numer == 0 | isnan(numer))
   r1 = NaN;
elseif a == 0 
   r1 = inf;
   
else
   r1 = numer/(2*a);
end
r2 = c/(a*r1);
tau = max(r1,r2); 
tau = min(1,tau);
if tau <= 0, 
   error('square root error in function deljust'); 
end
nx = tau*x;



