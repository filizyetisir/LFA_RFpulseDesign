


function sar2=prepare_sar(sar,spokes,system)
% function sar2=prepare_sar(sar,spokes,system)

sar2=sar;
nchannels=system.ncoils;

% local SAR matrices
file=fopen(sar2.lsarpath,'r');
tline=fgets(file); sar2.nvop=tline_to_param(tline,18,size(tline,2)-1);

if sar2.nvop<0 || sar2.nvop>1e4
    error('Error, # of VOP out of range in prepare_sar().');
end
tline=fgets(file); sar2.umax=tline_to_param(tline,51,size(tline,2)-1);
sar2.sarmats=zeros(sar2.nvop+1,nchannels,nchannels);
for n=1:sar2.nvop
    if n==1
        tline=fgets(file); tline=fgets(file);
    else
        tline=fgets(file); tline=fgets(file); tline=fgets(file);
    end
    for i=1:nchannels
        tmp=fscanf(file,'%e',2*nchannels);
        sar2.sarmats(n,i,:)=( tmp(1:2:2*nchannels-1)+1j*tmp(2:2:2*nchannels) ); %*spokes.sumsincsq/spokes.ntimes;
    end
    if norm( reshape(sar2.sarmats(n,:,:),nchannels,nchannels)-reshape(sar2.sarmats(n,:,:),nchannels,nchannels)' )>1e-7
        error( sprintf('Error, input local SAR matrix #%d is not hermitian in prepare_sar().',n) );
    end
end
fclose(file);

% global SAR matrix
file=fopen(sar2.gsarpath,'r');
for i=1:nchannels
    tmp=fscanf(file,'%e',2*nchannels);
    sar2.sarmats(sar2.nvop+1,i,:)=( tmp(1:2:2*nchannels-1)+1j*tmp(2:2:2*nchannels) ); %*spokes.sumsincsq/spokes.ntimes;
end
if norm( reshape(sar2.sarmats(sar2.nvop+1,:,:),nchannels,nchannels)-reshape(sar2.sarmats(sar2.nvop+1,:,:),nchannels,nchannels)' )>1e-7
    error( sprintf('Error, input global SAR matrix is not hermitian in prepare_sar().',n) );
end
fclose(file);



