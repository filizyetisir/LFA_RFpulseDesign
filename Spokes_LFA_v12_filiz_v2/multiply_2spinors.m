


function [a3 b3]=multiply_2spinors(a1,b1,a2,b2)
% function [a3 b3]=multiply_2spinors(a1,b1,a2,b2)
% perform spinor multiplication Q2*Q1 in this order

a3=a2.*a1 - conj(b2).*b1;
b3=b2.*a1 + conj(a2).*b1;

