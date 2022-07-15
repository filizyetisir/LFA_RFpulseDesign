

function [hess_an hess_numer]=test_hessian(fhandle,hesshandle,x)
% function [hess_an hess_numer]=test_hessian(fhandle,hesshandle,x)

n=size(x,1);

hess_an=hesshandle(x);

[f, df]=fhandle(x);

eps=1e-3;
hess_numer=zeros(n,n);

for i=1:n
    disp( sprintf('Computing Hessian row #%d out of %d ...',i,n) );
    for j=1:n
        if i==j
            xip=x; xip(i)=xip(i)+eps; [fip, dfip]=fhandle(xip);
            xim=x; xim(i)=xim(i)-eps; [fim, dfim]=fhandle(xim);
            hess_numer(i,j)=(fip+fim-2.0*f)/eps^2;
        else
            xipjp=x; xipjp(i)=xipjp(i)+eps; xipjp(j)=xipjp(j)+eps; [fipjp, dfipjp]=fhandle(xipjp);
            xipjm=x; xipjm(i)=xipjm(i)+eps; xipjm(j)=xipjm(j)-eps; [fipjm, dfipjm]=fhandle(xipjm);
            ximjp=x; ximjp(i)=ximjp(i)-eps; ximjp(j)=ximjp(j)+eps; [fimjp, dfimjm]=fhandle(ximjp);
            ximjm=x; ximjm(i)=ximjm(i)-eps; ximjm(j)=ximjm(j)-eps; [fimjm, dfimjm]=fhandle(ximjm);
            hess_numer(i,j)=(fipjp+fimjm-fipjm-fimjp)/(4.0*eps^2);
        end
    end
end

% save hess_numer.mat hess_numer;
% load hess_numer;

% figure; imagesc(hess_an); colormap(hot); colorbar; axis image; title('Hessian analytical'); % caxis([0 1e-7])
% figure; imagesc(hess_numer); colormap(hot); colorbar; axis image; title('Hessian numerical');  %caxis([0 1e-7])
% figure; imagesc(hess_an-hess_numer); colormap(hot); colorbar; axis image; title('Hessian analytical - Hessian numerical'); % caxis([0 1e-7])

figure; imagesc( [hess_an hess_numer hess_an-hess_numer ]); colormap(hot); colorbar; axis image;  title('An.   Num.   An.-Num.');

diff=abs( hess_an-hess_numer );
mean(diff(:))

diff=abs( hess_an-hess_numer ) ./ abs(hess_numer);
ind = (diff>1e9);
ind2 = isnan(diff);
diff(ind) = 0;
diff(ind2) = 0;
mean(diff(:))