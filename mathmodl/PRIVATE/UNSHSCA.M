function[xx] = unshsca(x,l,u,DS);
%UNSHSCA Unshift and unscale
%
%  xx = UNSHSCA(x,l,u,DS); vector x is shifted and scaled to yield xx.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1998/03/21 16:29:01 $

n = length(x);
arg1 = (l== -inf) & (u == inf);
arg2 = (l== -inf) & (u < inf); 
arg3 = (l> -inf) & (u == inf);
arg4 = (l > -inf) & (u < inf);
%
% UNSCALE
xx = full(DS*x);   % always full except in scalar case.
%
% UNSHIFT
xx(arg2) = xx(arg2) + u(arg2) -ones(size(arg2(arg2>0))); 
xx(arg3) = xx(arg3) + l(arg3);
xx(arg4) = xx(arg4) + l(arg4);


