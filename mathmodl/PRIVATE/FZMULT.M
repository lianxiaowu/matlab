function [W,L,U,pcol,P]  = fzmult(A,V,transpose,L,U,pcol,P);
%FZMULT Multiplication with fundamental nullspace basis.
%   W = FZMULT(A,V) Computes the product of a
%   a fundamental basis for the nullspace of full row-rank
%   matrix A, Z, times matrix V; W = Z*V.
%   Assume A is m-by-n with m<n, rank(A) = m.
%   Partition A = (A_1,A_2) where
%   A_1 is square. It is assumed that rank(A_1) = m.
%   The nullspace basis matrix Z
%   is not formed explicity. An implicit representation
%   is used based on the sparse LU-factorization of A_1.
%
%   W = FZMULT(A,V,TRANSPOSE), when TRANSPOSE is 'transpose'
%   computes the product of the transpose with V: W = Z'*V.
%
%   W = FZMULT(A,V,TRANSPOSE,L,U,pcol,P) uses the precomputed
%   sparse LU-factorization of matrix A_1 :  P * A_1(:,pcol) = LU.
%
%   [W,L,U,pcol,P]  = FZMULT(A,V) returns the sparse LU-factorization
%   of matrix A_1:  P* A_1(:,pcol) = LU

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1998/03/21 16:29:07 $

% Initialization
if nargin < 2
   error('function FZMULT expects at least 2 input arguments'); 
end
[m,n] = size(A); 
if m > n
   error('function FZMULT expects A to have more columnss than rows'); 
end
if nargin < 3 | isempty(transpose)
   transpose = '';     % default
end

[mm,p] = size(V);
switch transpose
case ''
   if mm ~= n-m
      error('Dimensions of A and V not compatible in function FZMULT'); 
   end
case 'transpose'
   if mm ~= n
      error('Dimensions of A and V not compatible in function FZMULT'); 
   end
otherwise
   error('Invalid string used for TRANSPOSE argument to FZMULT');
end

A1 = A(:,1:m);      % A1 is square 
A2 = A(:,m+1:n); 
if nargin  <  7
   pcol = colmmd(A1); 
   [L,U,P] = lu(A1(:,pcol));
end

switch transpose
case ''
   W = -A2*V; 
   WW = L\P*W; 
   W = U\WW; 
   WW(pcol,:) = W;
   W = [WW;V]; 
case 'transpose'
   V1 = V(1:m,:); 
   V2 = V(m+1:n,:);
   WW = -V1(pcol,:); 
   W = U'\WW; 
   WW = L'\W; 
   W = P'*WW; 
   WW = A2'*W;
   W = WW + V2; 
otherwise
   error('Invalid string used for TRANSPOSE argument to FZMULT');
end
W = full(W);






