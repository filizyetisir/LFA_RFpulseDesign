


function h=dualhess_STA(xrf,lambda,A,adj,spokes,system,sar,solver,do_mls,opt)
% function h=dualhess_STA(xrf,lambda,A,adj,spokes,system,sar,solver,do_mls)

[npix ncunk]=size(A);
nchannels=ncunk/spokes.nspokes;

% Hessian of the main function
rf=xrf(1:ncunk) + 1j*xrf(ncunk+1:end);
MT=A*rf;

if do_mls==1
    y=abs(MT) - abs(adj.target_STA);
    z=sqrt( abs(adj.target_STA)./abs(MT) ) .* exp(1j*angle(MT));
    w=y./abs(MT);
    
    z( abs(MT)==0 )=0;
    w( abs(MT)==0 )=0;
    
    % X=A'*diag(z);
    X=A'.*repmat(z.',[ncunk 1]);  % faster
    
    % Y=A'*diag(w)*A;
    Y=( A'.*repmat(w.',[ncunk 1]) )*A;  % faster
    
    h11=2*real(X)*real(X)' + 2*real(Y);
    h12=2*real(X)*imag(X)' - 2*imag(Y);
    h21=h12';
    h22=2*imag(X)*imag(X)' + 2*real(Y);
    
    h_OBJ=[ h11 h12;h21 h22 ];
else
    h11=2*A'*A;
    
    h_OBJ=[ real(h11)  -imag(h11);
            imag(h11)  real(h11) ];
end

% HESSIAN OF THE SAR AND POWER CONSTRAINTS
if solver==1  % Matlab solver
    LagMult=lambda.ineqnonlin;
else
    LagMult=lambda;
end

if(opt.sarcons)
    constant = spokes.sumsincsq/spokes.ntimes;
    % Hessian of the constraints (SAR)
    subh=zeros(system.ncoils,system.ncoils);
    nconst=sar.nvop+1;
    subh=2*squeeze( sum(constant*sar.sarmats .* repmat(LagMult(1:nconst),[1 nchannels nchannels]),1 ) );
    ch=kron( eye(spokes.nspokes),subh ).*(system.ncoils);
    h_SAR=[real(ch) -imag(ch)';imag(ch)' real(ch)];
end

% Hessian of the constraints (max. power)
if(opt.sarcons)
    tmp=LagMult(nconst+1:nconst+ncunk);
else
    tmp=LagMult(1:ncunk);
end
tmp=repmat(tmp,[2 1]); % real and imaginary
h_MAX_POW=2.0*diag(tmp);


% Added 8/16/2020
% Hessian of the constraints (avg power)
if(opt.sarcons)
    tmp=LagMult(nconst+ncunk+1:nconst+ncunk+nchannels);
else
    tmp=LagMult(ncunk+1: ncunk+nchannels);
end
tmp=repmat(tmp,[2*spokes.nspokes 1]); % multiple spokes and real and imaginary
h_AV_POW=2.0*diag(tmp);



if(opt.sarcons)
    h=h_OBJ + h_SAR + h_MAX_POW + h_AV_POW;
else
    h=h_OBJ + h_MAX_POW + h_AV_POW;
end





