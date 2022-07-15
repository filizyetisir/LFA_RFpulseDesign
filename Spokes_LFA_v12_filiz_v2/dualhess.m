

function h=dualhess(xrf,lambda,A,adj,spokes,sar)
% function h=dualhess(xrf,lambda,A,adj,spokes,sar)

nchannels=size(adj.b1maps,1);
ncunk=spokes.nspokes*nchannels;
% nconst=sar.nvop+1;

% Hessian of the main function
ch=2.0*A'*A;
h=[real(ch) -imag(ch);imag(ch) real(ch)];

% % Hessian of the constraints (SAR)
% subh=zeros(nchannels,nchannels);
% for n=1:nconst
%     subh(:,:)=subh(:,:) + 2.0*lambda.ineqnonlin(n)*reshape(sar.sarmats(n,:,:),nchannels,nchannels);
% end
% ch=kron( eye(spokes.nspokes),subh );
% h=h+[real(ch) -imag(ch)';imag(ch)' real(ch)];

% Hessian of the constraints (max. power)
% tmp=lambda.ineqnonlin(nconst+1:nconst+ncunk);
tmp=lambda.ineqnonlin(1:ncunk);
tmp=[tmp;tmp]; % real and imaginary
h=h+2.0*diag(tmp);

% % Hessian of the constraints (av. power)
% tmp=[];
% for i=1:spokes.nspokes
%     tmp=[tmp;lambda.ineqnonlin(nconst+ncunk+1:nconst+ncunk+nchannels)];
% end
% tmp=[tmp;tmp]; % real and imaginary
% h=h+2.0*diag(tmp);



