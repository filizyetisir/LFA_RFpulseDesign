

function rf=write_and_display_pulse(x,opt,spokes,system,title_str)
% function rf=write_and_display_pulse(x,opt,spokes,system,title_str)

nchannels=system.ncoils;

file=fopen( sprintf('%s/pulse.txt',opt.resdir) ,'w');
for j=1:spokes.nspokes
    fprintf(file,'spoke #%d mag.\tspoke #%d pha.\t',j,j);
end
for i=1:nchannels
    fprintf(file,'\n');
    for j=1:spokes.nspokes
		fprintf(file,'%e\t%e\t',abs(x(i+(j-1)*nchannels)),angle(x(i+(j-1)*nchannels))/pi*180.0);
    end
end

fprintf(file,'\ntime (s)\tgx (T/m)\tgy (T/m)\tgz (T/m)\t\tsx (T/m/s)\tsy (T/m/s)\tsz (T/m/s)\tkx (1/m)\tky (1/m)\tkz (1/m)');
for j=1:nchannels
    fprintf(file,'\trf_%d_mag (V)\trf_%d_pha (deg.)',j,j);
end
time=[1:spokes.ntimes]'*system.deltat;
rf=zeros(spokes.ntimes,nchannels);
for i=1:spokes.ntimes
    fprintf(file,'\n');
    if i==1
        sr=[0 0 0];
    else
        sr=(spokes.g(i,:)-spokes.g(i-1,:))/system.deltat;
    end
	fprintf(file,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e',time(i),spokes.g(i,:),sr,spokes.k(i,:));
    for j=1:nchannels
        rf(i,j)=get_rf(x,i,j,nchannels,spokes);
        fprintf(file,'\t%e\t%e',abs(rf(i,j)),angle(rf(i,j))/pi*180.0);
    end
end
fclose(file);

if strcmpi(opt.display,'on')
    figure;
    subplot(2,1,1); plot(time,abs(rf)); xlabel('time [s]'); title(sprintf('RF mag. (%s) [V]',title_str))
    subplot(2,1,2); plot(time,angle(rf)/pi*180.0); xlabel('time (s)'); title(sprintf('RF pha. (%s) [V]',title_str))
end
    

% write ini file
options.Oversampling=1;
options.NominalFlipAngle=opt.tfa;
% options.NOMFLIPANGLE=opt.tfa;
%options.NOMFLIPANGLE=15;  % although the actual flip-angle is defined by opt.tfa, it is useful to "pretend" that the 
                           % it is always 15 degrees, which is the
                           % flip-angle in the default FLASH_ptx protocol.
                           % This avoids having to manually change the
                           % protocol flip-angle so that it matches the
                           % designed flip-angle (if there is a mismatch between the two, then
                           % the pulse is scaled by their ratio in the
                           % sequence)

options.RFPULSE_ID=0;
options.VERBOSE=0;
options.GRADRASTERTIME=system.deltat*1e6;  % us
grad=[spokes.g(:,2) -spokes.g(:,1) -spokes.g(:,3)]*1e3;  % swap gx and gy and convert T/m -> mT/m (the delayed version of the gradienst is sent to the scanner)

% save ini file
grad_over_samp_fac=1e-5/system.deltat;
if ceil(grad_over_samp_fac)~=grad_over_samp_fac
    warning('The oversampling factor is not integer, this will cause the ini pulse to be unplayable on the scanner.');
else
    options.GRADRASTERTIME=system.deltat*grad_over_samp_fac*1e6;  % 10 us
    options.FACTOROVERSAMPLE=grad_over_samp_fac;
    grad=grad(1:grad_over_samp_fac:end,:);
end

% save_pTXRFPulse_toINI(grad,reshape(rf,system.ncoils*opt.ntimes,1),[],options);

options.type= 'SBBExcitationPtx';
options.fileName='pTXArbitrary';
options.PulseName= 'pTx';
options.Comment = 'pTx pulse';
create_pTXRFPulse2( grad.',rf', options);  % save the conjugate of the RF pulse

   
 
    
function rf=get_rf(x,time,channel,nchannels,spokes)
% function rf=get_rf(time,channel)

rf=0.0;
if time>0 && time<=spokes.ntimes && spokes.time_to_spoke(time)>0 && channel>0 && channel<=nchannels
    rf=x(channel + (spokes.time_to_spoke(time)-1)*nchannels)*spokes.sinc(spokes.time_to_sinc_time(time));
end







    
    
    
    
    
    
    
    
    
    
    
    
