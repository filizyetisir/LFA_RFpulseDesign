function r = rect (x)
% RECT   The rectangular function is defined to be 1
%        on [-0.5,0.5] and 0 elsewhere
%
% m-file function [rect.m] devised for practice midterm 2 in cs 205/fall 96
% the file is available from the class home page
%     http://www.stewart.cs.sdsu.edu/cs205/
% and from the individual student's rohan account for this course
% so that you can practice with it to check your understanding of
% communications between MATLAB m-file scripts and m-file functions
%
% initialize output matrix to be all zeros, same size as input
r = zeros(size(x)); 
%
set1 = find(abs(x) <= 0.5); % identify component of input vector in [-0.5,0.5]
r(set1) = ones(size(set1)); % for those component in [-0.5,0.5], set value 1
