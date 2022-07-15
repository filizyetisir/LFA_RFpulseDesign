
function spokes2=prepare_gradient(spokes,system)
% function spokes2=prepare_gradient(spokes,system)

spokes2=spokes;
spokes2.kextent=spokes.nsinczc/spokes.sthick;

gamma=42.576*1e6;

spokes2.sdir=spokes2.sdir/norm(spokes2.sdir);
sdir2=spokes2.sdir;

p1=spokes2.spcoords(1,:)'*(-1) + sdir2*( spokes2.kextent/2.0 );
direc=p1; kdist=norm( direc );
grad=translate_kspace(kdist,direc,spokes,system);
t2s=-1*ones( size(grad,1),1 );
t2st=-1*ones( size(grad,1),1 );

for sp=1:spokes2.nspokes
    % spoke move
    p2=spokes2.spcoords(sp,:)'*(-1) - sdir2*( spokes2.kextent/2.0 );
	direc=p2-p1; kdist=norm( direc );
    tmp=translate_kspace(kdist,direc,spokes,system);
    grad=[ grad;tmp ];
	if sp==1
        spokes2.ntimesperspoke=size( tmp,1 ); % # of time points per spoke
    end
    t2s=[ t2s;sp*ones( size(tmp,1),1 ) ];
    t2st=[ t2st;[1:size(tmp,1)]' ];
    
    % stitch move
    p1=p2;
    if sp<spokes2.nspokes
        p2=spokes2.spcoords(sp+1,:)'*(-1) - sdir2*( spokes2.kextent/2.0 );
        direc=p2-p1; kdist=norm( direc );
        if kdist>1e-10
            tmp=translate_kspace(kdist,direc,spokes,system);
        else
            tmp=[0 0 0];
        end
        grad=[ grad;tmp ];
        p1=p2;
        sdir2=-1*sdir2;
    end
    t2s=[ t2s;-1*ones( size(tmp,1),1 ) ];
    t2st=[ t2st;-1*ones( size(tmp,1),1 ) ];
end

% compute trajectory
spokes2.ntimes=size(grad,1);
spokes2.length=(spokes2.ntimes)*(system.deltat);
spokes2.g=zeros(spokes2.ntimes,3); spokes2.k=zeros(spokes2.ntimes,3);
spokes2.time_to_spoke=zeros(spokes2.ntimes,1);
spokes2.time_to_sinc_time=zeros(spokes2.ntimes,1);
for i=1:spokes2.ntimes
    if t2s(spokes2.ntimes-i+1)>=0
        spokes2.time_to_spoke(i,1)=spokes2.nspokes-t2s(spokes2.ntimes-i+1)+1;
        spokes2.time_to_sinc_time(i,1)=spokes2.ntimesperspoke-t2st(spokes2.ntimes-i+1)+1;
    else
        spokes2.time_to_spoke(i,1)=-1;
        spokes2.time_to_sinc_time(i,1)=-1;
    end
    spokes2.g(i,:)=grad(spokes2.ntimes-i+1,:);
end

for i=2:spokes2.ntimes
    spokes2.k(i,:)=spokes2.k(i-1,:) + ( grad(i,:)+grad(i-1,:) )/2.0*(system.deltat)*gamma;
end
spokes2.k(:,:)=-spokes2.k(end:-1:1,:);

% check that all spokes have the same # of samples
spokes2.spokes_time_offsets=zeros(spokes2.nspokes,1);
j1=1; j2=1; j=0;
for i=1:spokes2.nspokes
    spokes2.spokes_time_offsets(i,1)=j;
    j2=j1;
    while spokes2.time_to_spoke(j1,1)>=0
        j1=j1+1; j=j+1;
    end
    if (j1-j2)~=spokes2.ntimesperspoke
        error('Error, not all spokes have the same # of samples in prepare_gradient().');
    end
    
    while j1<=spokes2.ntimes && spokes2.time_to_spoke(j1,1)<0
        j1=j1+1; j=j+1;
    end
end

spokes2.sinc_time_and_spoke_to_time=zeros(spokes2.ntimesperspoke,spokes2.nspokes);
for sp=1:spokes2.nspokes
    for sinctime=1:spokes2.ntimesperspoke
        spokes2.sinc_time_and_spoke_to_time(sinctime,sp)=get_time_from_spoke_and_sinctime(spokes2,sp,sinctime);
    end
end



function grad=translate_kspace(ktrans,direc,spokes,system)

gamma=42.576*1e6;

gmax=1*(system.gmax);  % downgrade gradient performance here
smax=1*(system.smax);

trise=gmax/smax;
ntri=floor(trise/system.deltat) + 1;
ktri=0.5*gmax*ntri*(system.deltat)*gamma;
if ktri<ktrans/2
    % flat top required
    kftop=ktrans - 2.0*ktri;
    nftop=floor(kftop/gamma/gmax/system.deltat);
    n=ntri*2+nftop+1;
    g=zeros(n,1);
    for i=1:ntri
        g(i+1)=i/ntri*gmax;
    end
    for i=1:nftop
        g(i+1+ntri)=gmax;
    end
    for i=1:ntri
        g(i+1+ntri+nftop)=gmax - i/ntri*gmax;
    end
    
else
    % no flat top required
    tramp=sqrt( ktrans*ntri*system.deltat/gamma/gmax );
    nramp=floor(tramp/system.deltat)+1;
    gpeak=gmax*nramp/ntri;
    n=2*nramp+1;
    g=zeros(n,1);
    for i=1:nramp
        g(i+1)=i/nramp*gpeak;
    end
    for i=1:nramp
        g(i+1+nramp)=gpeak - i/nramp*gpeak;
    end
end

% adapt gradient amplitude
k=zeros(n,1);
for i=2:n
    k(i)=k(i-1) + ( g(i)+g(i-1) )/2.0*system.deltat*gamma;
end
for i=1:n
    g(i)=g(i)*ktrans/abs(k(n));
end

% rotate gradient waveform
direc=direc/norm( direc );
for i=1:n
    grad(i,1)=g(i)*direc(1);
    grad(i,2)=g(i)*direc(2);
    grad(i,3)=g(i)*direc(3);
end

clear g; clear k;




