function [f,F,V]=Companion(B,V,N,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		Luca Benati
%		European Central Bank
%		Monetary Policy Strategy Division
%		January 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f,F,V]=Companion(B,V,N,p)
% This program computes the companion form of a VAR(p).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              Input of the program:
%
% B = the key VAR matrix: B = [B(0) B(1) B(2) ... B(p)]
% V = the VAR's variance-covariance matrix
% N = the number of variables in the VAR
% p = the order of the VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              Output of the program:
% f      = the intercept vector of the companion form
% F      = the companion matrix
% V      = the companion form's variance-covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B0=B(:,1);
f=[B0; zeros(N*(p-1),1)];
B=B(:,2:size(B,2));
F=[B; eye(N*(p-1)) zeros(N*(p-1),N)];
V=[V zeros(N,N*(p-1)); zeros(N*(p-1),N) zeros(N*(p-1),N*(p-1))];












