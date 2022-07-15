


function [a4 b4]=multiply_3spinors(a1,b1,a2,b2,a3,b3)
% function [a4 b4]=multiply_3spinors(a1,b1,a2,b2,a3,b3)
% perform the spinor product Q3*Q2*Q1 in this order

tmp1=a2.*a1 - conj(b2).*b1;
tmp2=b2.*a1 + conj(a2).*b1;

a4=a3.*tmp1 - conj(b3).*tmp2;
b4=b3.*tmp1 + conj(a3).*tmp2;









