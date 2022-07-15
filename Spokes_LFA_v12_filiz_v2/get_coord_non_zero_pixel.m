

function [i j k x y z]=get_coord_non_zero_pixel(opt,n)
% function [i j k x y z]=get_coord_non_zero_pixel(opt,n)

i=opt.nonzeropixels(n,1);
j=opt.nonzeropixels(n,2);
k=opt.nonzeropixels(n,3);
x=opt.nonzeropixels(n,4);
y=opt.nonzeropixels(n,5);
z=opt.nonzeropixels(n,6);

