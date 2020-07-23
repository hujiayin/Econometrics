%%% A function to generate Hermite using a loop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [temp]=Hermite(k)
% This function computes the Hermite Polynomial Recursively
% Store this function in M-file Hermite.m
% Note Hermite(1)=1

syms z

H{1}=sym('1');

for n=2:k
    H{n}=simplify(z*H{n-1}-diff(H{n-1},z));
end

temp=H{k};
