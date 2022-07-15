
function img=conv_3Dstack_2_2Dimg(stack,nfigsx,nfigsy)
% function img=conv_3Dstack_2_2Dimg(stack,nfigsx,nfigsy)

[nx ny nz]=size(stack);

img=[];
n=1;
for i=1:nfigsx
    row=[];
    for j=1:nfigsy
        if n<=nz
            row=[ row stack(:,:,n) ];
        else
            row=[ row zeros(nx,ny) ];
        end
        n=n+1;
    end
    img=[ img;row ];
end



