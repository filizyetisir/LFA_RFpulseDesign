

function [c ceq dc dceq]=fn(xrf,adj,spokes,sar,system,opt)
% function [c ceq dc dceq]=fn(xrf,adj,spokes,sar,system)

% no equality constraints
ceq=[];
dceq=[];

% inequality constraints (SAR)
nchannels=system.ncoils;
if(opt.sarcons)
    nconst=sar.nvop+1;
end
%nconst=sar.nvop+1;
ncunk=size(xrf,1)/2;


constant = spokes.sumsincsq/spokes.ntimes;
    
    
rf=reshape( xrf(1:ncunk) + 1j*xrf(ncunk+1:2*ncunk),1,ncunk );

c=[];
dc=[];

if(opt.sarcons)
    c=zeros(nconst,1);
    for i=1:nconst
        if i==nconst
            c(i)=constant*real( rf*(kron(eye(spokes.nspokes),reshape(sar.sarmats(i,:,:),nchannels,nchannels)))*rf' ) - sar.gsarmax;
        else
            c(i)=constant*real( rf*(kron(eye(spokes.nspokes),reshape(sar.sarmats(i,:,:),nchannels,nchannels)))*rf' ) - sar.lsarmax;
        end
    end
    
    % inequality constraints gradient (SAR)
    dc=zeros(2*ncunk,nconst);
    for i=1:nconst
        tmp=2.0*constant*(kron(eye(spokes.nspokes),reshape(sar.sarmats(i,:,:),nchannels,nchannels)))*rf';
        dc(1:ncunk,i)=real(tmp);
        dc(ncunk+1:2*ncunk,i)=-imag(tmp);
    end
end

% inequality constraints (max. power)
% xsqmax=sar.ppmax*8.0*50.0;  % P=U^2/(8*Z0)
xsqmax=system.umax^2;
c=[ c;abs(rf.').^2-xsqmax ];
dc=[ dc [2*real(diag(rf(1,:)));2*imag(diag(rf(1,:)))] ];


% % inequality constraints (av. power)
rf2=reshape(rf.',nchannels,spokes.nspokes).';
xsumsqmax=sar.ppavmax*8.0*50.0*(spokes.ntimes/spokes.sumsincsq);
c=[ c;sum(abs(rf2).^2,1)'-xsumsqmax ];
tmp2=[];
for i=1:spokes.nspokes
    tmp2=[ tmp2;diag(rf2(i,:)) ]; 
end
dc=[dc [2*real(tmp2);2*imag(tmp2)]];

