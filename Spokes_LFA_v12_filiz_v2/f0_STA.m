


function [f df MT]=f0_STA(xrf,A,adj,do_mls)
% function [f df MT]=f0_STA(xrf,A,adj,do_mls)

% function
ncunk=size(xrf,1)/2;
rf=xrf(1:ncunk,1) + 1j*xrf(ncunk+1:2*ncunk,1);
MT=A*rf;

if do_mls==1
    y=abs(MT) - abs(adj.target_STA);  % better, avoids the phase curl problem
else
    y=MT - adj.target_STA;
end
f=(y'*y);

% gradient
if nargout>1        
    if do_mls==1
        cdf=2.0*A'*( y.*exp(1j*angle(MT)) );
    else
        cdf=2.0*A'*y;
    end    
    df(1:ncunk,1)=real(cdf);
    df(ncunk+1:2*ncunk,1)=imag(cdf);
end






