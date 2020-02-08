function[val,g] = fquad(x,c,H,mtxmpy,data,D)
%FQUAD Evaluate quadratic function.
%	val = FQUAD(x,c,H,mtxmpy,data,D) evaluates the quadratic
%   function val = c'*x + .5*x'*D*MTX*D*x, where
%   D is a diagonal matrix and MTX is defined by the 
%   matrix multiply routine 'mtxmpy' and 'data'.
%
%   [val,g] = FQUAD(x,c,H,mtxmpy,data,D) also evaluates the 
%   gradient: g = D*MTX*D*x + c.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1998/03/21 16:29:06 $

if nargin < 4, 
   error('fquad requires at least 4 input parameters'), 
end
n = length(x);
if nargin < 5, 
   data = []; 
end
if nargin < 6, 
   D = speye(n); 
end
w = full(D*x);                  % w always full unless both scalar
ww = feval(mtxmpy,w,H,data);
w = full(D*ww);                 % w always full unless both scalar
g = w + c; 
val = x'*((0.5*w) + c);


