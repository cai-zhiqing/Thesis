function InvX=MyInverse(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luca Benati
% June 24, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function InvX=MyInverse(X)
% This function computes an inverse for a (KxK) matrix X. The key reason for
% using this function, rather than the MATLAB built-in function inv.m,
% is that in the 2014 version of MATLAB, inv.m performs very poorly,
% producing unreliable results. 
%                                             Input of the program is:
% X    = a square matrix
%                                             Output of the program is:
% InvX = its inverse
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[K,K2]=size(X);
%
if K~=K2
    disp('Matrix is not square: Exit')
    InvX=NaN;
    return
end
% The inverse:
InvX=X\eye(K);
