


function [mxf, myf, mzf, beta2, mxy2D_timely]=bloch_simulation(rf,adj,spokes,system,freq_offset,M0)
% function [mxf, myf, mzf, beta2, mxytimely]=bloch_simulation(rf,adj,spokes,system,freq_offset,M0)

% Simulate the pulse for the whole data set 
adj.roi = adj.roiWhole;
adj.b0map = adj.b0mapWhole;
adj.b1maps = adj.b1mapsWhole;
adj.ind=find(adj.roi>0);
adj.nnonzeropixels=size(adj.ind,1);
[indx, indy, indz]=ind2sub( [adj.nx adj.ny adj.nz],adj.ind );
[tmpx, tmpy, tmpz]=ndgrid(adj.xs,adj.ys,adj.zs);
adj.nonzeropixels=[indx indy indz tmpx(adj.ind) tmpy(adj.ind) tmpz(adj.ind)];


% initial magnetization state
mxy0=M0(1) + 1j*M0(2);
mz0=M0(3);  

gamma=42.576*1e6;
nchannels=system.ncoils;

fprintf('Bloch simulation (0%%) ...');
mxf=zeros(adj.nx,adj.ny,adj.nz);
myf=zeros(adj.nx,adj.ny,adj.nz);
mzf=zeros(adj.nx,adj.ny,adj.nz);
beta2=zeros(adj.nx,adj.ny,adj.nz);
q=zeros(adj.nnonzeropixels,4); q(:,1)=1; q(:,4)=1;
b1maps=zeros(nchannels,adj.nnonzeropixels);
ind=sub2ind( [adj.nx adj.ny adj.nz], adj.nonzeropixels(:,1),adj.nonzeropixels(:,2),adj.nonzeropixels(:,3) );
for chn=1:nchannels
    b1map=reshape( adj.b1maps(chn,:,:,:),adj.nx,adj.ny,adj.nz );
    b1maps(chn,:)=b1map(ind);
end
coords=adj.nonzeropixels(:,4:6);

% ntimes=ceil( 10e-3 / system.deltat );  % force 10 ms simulation
ntimes=spokes.ntimes;

percent=10;
% for time=1:spokes.ntimes
for time=1:ntimes
    
    % progress
    if mod(time,floor(ntimes/10))==0
        if percent==10
            fprintf('\b\b\b\b\b\b\b\b(%d%%) ...',percent);
        elseif percent==100
            fprintf('\b\b\b\b\b\b\b\b\b(%d%%) ...\n',percent);
        else
            fprintf('\b\b\b\b\b\b\b\b\b(%d%%) ...',percent);
        end
        percent=percent+10;
    end

    % Cayley-Klein parameters
    if time<=spokes.ntimes
        b1s=( rf(time,:)*b1maps ).';
        grad=coords*spokes.g(time,:)' - ( adj.b0map(ind)+freq_offset )/gamma;
    else
        b1s=zeros(size(b1s));
        grad=-( adj.b0map(ind)+freq_offset )/gamma;
    end    
    norm=sqrt( abs(b1s).^2+grad.^2 );
    phi=-2*pi*gamma*system.deltat*norm;
    if sum(abs(phi))<=0.0
        continue;
    end
    n=[real(b1s) imag(b1s) grad];
    n=n./repmat(norm,1,3);
    a=cos(phi/2.0) - 1j*(n(:,3)).*sin(phi/2.0);
    b=sin(phi/2.0).*n(:,2) - 1j*sin(phi/2.0).*n(:,1);    
    ind2=find(norm==0 | isnan(norm));
    a(ind2)=1.0;
    b(ind2)=0.0;
    
    % update rotation matrix
    q1=q(:,1); q2=q(:,2); q3=q(:,3); q4=q(:,4);
    q(:,1)=a.*q1 - conj(b).*q3;
    q(:,2)=a.*q2 - conj(b).*q4;
    q(:,3)=b.*q1 + conj(a).*q3;
    q(:,4)=b.*q2 + conj(a).*q4;
    
    atemp=q(:,1); btemp=q(:,3);
    temp = zeros(adj.nx,adj.ny,adj.nz);
    temp(ind) =(conj(atemp).*conj(atemp))*mxy0 - (btemp.*btemp).*conj(mxy0) + (2.0*conj(atemp).*btemp).*mz0;
    mxy2D_timely(:,:,:,time) = temp;
    
end

% final magnetization
a=q(:,1); b=q(:,3);
mxy=(conj(a).*conj(a))*mxy0 - (b.*b).*conj(mxy0) + (2.0*conj(a).*b).*mz0;
mz=(-conj(a).*conj(b))*mxy0 - (a.*b).*conj(mxy0) + (a.*conj(a)-b.*conj(b))*mz0;

mxf(ind)=real(mxy);
myf(ind)=imag(mxy);
mzf(ind)=real(mz);
beta2(ind)=b.^2;







