
function spokes=prepare_gradient_spokes(spokes,system,opt)
% function spokes=prepare_gradient(spokes,system)

gamma=42.576*1e6;

Gss=(system.gmax)*(spokes.SS_grad_down);
spokes.kextent=spokes.nsinczc/spokes.sthick + gamma*Gss^2/system.smax;

spokes.sdir=spokes.sdir/norm(spokes.sdir);
sdir2=spokes.sdir;

p1=spokes.spcoords(1,:)'*(-1) + sdir2*( spokes.kextent/2.0 );
direc=p1; kdist=norm( direc );
grad=translate_kspace(kdist,direc,spokes,system,0);
t2s=-1*ones( size(grad,1),1 );
t2st=-1*ones( size(grad,1),1 );
mask_rf=zeros( size(grad,1),1 );

for sp=1:spokes.nspokes
    
    % spoke move
    p2=spokes.spcoords(sp,:)'*(-1) - sdir2*( spokes.kextent/2.0 );
	direc=p2-p1; kdist=norm( direc );
    [grad_ mask_rf_]=translate_kspace(kdist,direc,spokes,system,1);
    
    % create RF (possibly using VERSE)
    [spokes grad_ mask_rf_]=create_rf(spokes,system,grad_,mask_rf_, opt);

    % update spokes structures
    grad=[ grad;grad_ ];
    mask_rf=[ mask_rf;mask_rf_ ];
    t2s=[ t2s;sp*ones( size(grad_,1),1 ) ];
    t2st=[ t2st;[1:size(grad_,1)]' ];
	if sp==1
        spokes.ntimesperspoke=size( grad_,1 ); % # of time points per spoke
    end

    % stitch move
    p1=p2;
    if sp<spokes.nspokes
        
        % p2=spokes.spcoords(sp+1,:)'*(-1) - sdir2*( spokes.kextent/2.0 );
        p2=spokes.spcoords(sp+1,:)'*(-1) + sdir2*( spokes.kextent/2.0 );
        
        direc=p2-p1; kdist=norm( direc );
        if kdist>1e-10
            tmp=translate_kspace(kdist,direc,spokes,system,0);
        else
            tmp=[0 0 0];
        end
        grad=[ grad;tmp ];
        p1=p2;
        
        % sdir2=-1*sdir2;
        
        t2s=[ t2s;-1*ones( size(tmp,1),1 ) ];
        t2st=[ t2st;-1*ones( size(tmp,1),1 ) ];
        mask_rf=[ mask_rf;zeros( size(tmp,1),1 ) ];
        
    elseif(sp == spokes.nspokes)
        if(opt.gsym == 1)
        % Final stictch 12/16/16
        % This is to make the whole trajectory symmetrically balanced
        % as in Kawins paper rather than having a rewinder right after
        p2=[0 0 0]';
        direc=p2-p1; kdist=norm( direc );
        if kdist>1e-10
            tmp=translate_kspace(kdist,direc,spokes,system,0);
        else
            tmp=[0 0 0];
        end
        grad=[ grad;tmp ];
        p1=p2;
        
        t2s=[ t2s;-1*ones( size(tmp,1),1 ) ];
        t2st=[ t2st;-1*ones( size(tmp,1),1 ) ];
        mask_rf=[ mask_rf;zeros( size(tmp,1),1 ) ];
        end
    
    end
end


spokes.mask_rf=mask_rf(end:-1:1,1);

% compute trajectory
spokes.ntimes=size(grad,1);
spokes.length=(spokes.ntimes)*(system.deltat);
spokes.g=zeros(spokes.ntimes,3); 
spokes.time_to_spoke=zeros(spokes.ntimes,1);
spokes.time_to_sinc_time=zeros(spokes.ntimes,1);
for i=1:spokes.ntimes
    if t2s(spokes.ntimes-i+1)>=0
        spokes.time_to_spoke(i,1)=spokes.nspokes-t2s(spokes.ntimes-i+1)+1;
        spokes.time_to_sinc_time(i,1)=spokes.ntimesperspoke-t2st(spokes.ntimes-i+1)+1;
    else
        spokes.time_to_spoke(i,1)=-1;
        spokes.time_to_sinc_time(i,1)=-1;
    end
    spokes.g(i,:)=grad(spokes.ntimes-i+1,:);
end



ind = 0;
%Eliminate rewinding gradient for 180, only works for 1 spoke for now
% if(opt.refoc_pulse == 1)
% 
%     for i = size(spokes.g,1):-1:1
%         if(spokes.time_to_sinc_time(i) > 0)
%             break;
%         end
%         ind = ind+1;
%     end
%     % Remove rewinder from the end
%     ncut = size(spokes.g,1)-ind+1;
%     gnew = spokes.g(1:ncut,:);
%     ntimesnew = size(gnew,1);
%     time_to_spoke_new = spokes.time_to_spoke(1:ncut);
%     time_to_sinc_time_new = spokes.time_to_sinc_time(1:ncut);
%     mask_rf_new = spokes.mask_rf(1:ncut);
%     ind = 0;
%    %  Add rewinder to the beginning
% %     rewinder = spokes.g(end-ind+1:end,:);
% %     gnew = [rewinder; spokes.g];
% %     ntimesnew = size(gnew,1);
% %     L = size(rewinder,1);
% %     time_to_spoke_new = [-1*ones(L,1); spokes.time_to_spoke];
% %     time_to_sinc_time_new = [-1*ones(L,1); spokes.time_to_sinc_time];
% %     mask_rf_new = [zeros(L,1); spokes.mask_rf];
% %     
%     spokes.g = gnew;
%     spokes.ntimes = ntimesnew;
%     spokes.time_to_spoke = time_to_spoke_new;
%     spokes.time_to_sinc_time = time_to_sinc_time_new;
%     spokes.mask_rf = mask_rf_new;
%     
%     gradnew = spokes.g(end:-1:1,:);
%     grad = gradnew;
%     
% end


spokes.k=zeros(spokes.ntimes,3);
for i=2:spokes.ntimes
    spokes.k(i,:)=spokes.k(i-1,:) + ( grad(i,:)+grad(i-1,:) )/2.0*(system.deltat)*gamma;
end
spokes.k(:,:)=-spokes.k(end:-1:1,:);


if(opt.gsym == 1)
    % check that all spokes have the same # of samples
    spokes.spokes_time_offsets=zeros(spokes.nspokes,1);
    j1=1; j2=1; j=0;
    for i=1:spokes.nspokes
        while j1<=spokes.ntimes && spokes.time_to_spoke(j1,1)<0
            j1=j1+1; j=j+1;
        end
        spokes.spokes_time_offsets(i,1)=j1;
        
        j2=j1;
        while spokes.time_to_spoke(j1,1)>=0
            j1=j1+1; j=j+1;
        end
        
        if (j1-j2)~=spokes.ntimesperspoke
            error('Error, not all spokes have the same # of samples in prepare_gradient().');
        end
    end
    
else
    % Use when the gradient is not symmetrically balanced
    j1=ind+1; j2=ind+1; j=0;
    for i=1:spokes.nspokes
        spokes.spokes_time_offsets(i,1)=j+ind;
        j2=j1;
        while spokes.time_to_spoke(j1,1)>=0
            j1=j1+1; j=j+1;
        end
        if (j1-j2)~=spokes.ntimesperspoke
            error('Error, not all spokes have the same # of samples in prepare_gradient().');
        end
    
        while j1<=spokes.ntimes && spokes.time_to_spoke(j1,1)<0
            j1=j1+1; j=j+1;
        end
    end
    
end

spokes.sinc_time_and_spoke_to_time=zeros(spokes.ntimesperspoke,spokes.nspokes);
for sp=1:spokes.nspokes
    for sinctime=1:spokes.ntimesperspoke
        spokes.sinc_time_and_spoke_to_time(sinctime,sp)=get_time_from_spoke_and_sinctime(spokes,sp,sinctime);
    end
end



function [grad mask_rf]=translate_kspace(ktrans,direc,spokes,system,downgradegrad)

gamma=42.576*1e6;

if downgradegrad==1
    gmax=system.gmax*spokes.SS_grad_down;  % downgrade gradient performance here for variable spoke length
else
    gmax=system.gmax;
end
smax=system.smax;

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
    
    mask_rf=[ zeros(ntri,1);ones(nftop+1,1);zeros(ntri,1) ];  % only play RF during flat top
    
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
    
    mask_rf=zeros(n,1);
    
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



function [spokes grad mask_rf]=create_rf(spokes,system,grad,mask_rf,opt)

gamma=42.576*1e6;

g=grad(:,3);  
spokes.Gss=max(abs(g));  % slice selection plateau
omega=gamma*(spokes.Gss)*(spokes.sthick);

N=size(grad,1);

times=([1:N]'-0.5-N/2)*system.deltat;

ind_rf=find(mask_rf(1:N)>0);
Nrf=ind_rf(end)-ind_rf(1);
Trf=Nrf*(system.deltat);

x=omega*times;        

f1=sinc(x);
f2=rect( times/Trf );    
f3=0.54 + 0.46*cos( 2*pi*times/Trf );  % hamming window
% f3=0.50 + 0.50*cos( 2*pi*times/Trf );  % hann window
        
rf=f1.*f2.*f3;

if spokes.verse_factor>0  % VERSE
    % fprintf('************ VERSEing THE SPOKE SUB-PULSE ************\n');
    
    if(opt.use_fixed_subpulse == 1)
        if exist('verse_subpulse.mat','file')
            load verse_subpulse;
        else
            Gz=abs([0;abs(g);0]);
            rf=[0;rf;0];
            
            Gmax=( system.gmax )*0.9;
            Smax=( system.smax )*0.9;
            [rfv gv]=mintverse(rf,Gz,system.deltat,spokes.verse_factor,Gmax,Smax,system.deltat);
            
            figure;
            subplot(2,1,1); plot(rfv);
            subplot(2,1,2); plot(gv);
            
            save verse_subpulse.mat rfv gv;
        end
    else
        Gz=abs([0;abs(g);0]);
        rf=[0;rf;0];
            
        Gmax=( system.gmax )*0.9;
        Smax=( system.smax )*0.9;
        [rfv gv]=mintverse(rf,Gz,system.deltat,spokes.verse_factor,Gmax,Smax,system.deltat);
        save verse_subpulse.mat rfv gv;
    end
    
    spokes.sinc=rfv;
    
    Nv=size(rfv,1);    
    gv=gv * sign(g(ceil(N/2)));
    grad=[ zeros(Nv,2) gv ];  
    mask_rf=abs(rfv)>1e-5;
else
    spokes.sinc=rf;
end
















