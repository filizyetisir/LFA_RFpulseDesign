

function adj2=prepare_optimization(opt,spokes,system,adj)
% function adj2=prepare_optimization(opt,spokes,system,adj)

adj2=adj;


% Decide the roi,b1map and b0map to use
if(opt.data == 1)
    adj2.roi = adj.roiWhole;
    adj2.b1maps = adj.b1mapsWhole;
    adj2.b0map = adj.b0mapWhole;
elseif(opt.data == 2)
    adj2.roi = adj.roiInner;
    adj2.b1maps = adj.b1mapsInner;
    adj2.b0map = adj.b0mapInner;
elseif(opt.data == 3)
    adj2.roi = adj.roiSpecial;
    adj2.b1maps = adj.b1mapsSpecial;
    adj2.b0map = adj.b0mapSpecial;
else
    error('Data type unrecognized')
end

% % non-zero pixels
% if strcmpi(opt.display,'on')
%     figure; imagesc(adj2.roi); axis image off; title('Optimization ROI');
% end

adj2.ind=find(adj2.roi>0);
adj2.nnonzeropixels=size(adj2.ind,1);
[indx indy indz]=ind2sub( [adj2.nx adj2.ny adj2.nz],adj2.ind );
[tmpx tmpy tmpz]=ndgrid(adj2.xs,adj2.ys,adj2.zs);
adj2.nonzeropixels=[indx indy indz tmpx(adj2.ind) tmpy(adj2.ind) tmpz(adj2.ind)];

% optimization targets -- STA
if opt.do_mls==1
    adj2.targetImg_STA=repmat(adj2.roi .* sin( 30/180.0*pi ), [length(opt.freqs),1,1]);  % STA design is always 30 degrees that is scaled afterward
    adj2.target_STA=repmat(adj2.roi(adj2.ind) .* sin( 30/180.0*pi ), [length(opt.freqs),1]);
else
    tmp = []; tmp2 = [];
    if opt.refoc_pulse == 1
        for i = 1:opt.nfreqs
            ttt=adj2.roi .* sin( 30/180.0*pi ) .* exp( 1j*2*opt.tphase(:,:,i) ).*(-1);
            tmp = [tmp;ttt];
            tmp2 = [tmp2; ttt(adj2.ind)]; 
        end
        adj2.targetImg_STA = tmp;
        adj2.target_STA = tmp2;
    else
        tmp = []; tmp2 = [];
        for i = 1:opt.nfreqs
            ttt = adj2.roi .* sin( 30/180.0*pi ) .* exp( 1j*opt.tphase(:,:,i) );
            tmp = [ tmp; ttt];
            tmp2 = [tmp2; ttt( adj2.ind)];
        end
        adj2.targetImg_STA = tmp;
        adj2.target_STA = tmp2;
    end
end
%adj2.target_STA=adj2.targetImg_STA( adj2.ind );


% optimzation targets -- LFA
tfa_rad=opt.tfa/180.0*pi;
if opt.refoc_pulse==0  % excitation pulse   
    if opt.do_mls==1  % MLS
        adj2.targetImg_LFA=repmat(adj2.roi .* cos(tfa_rad),[length(opt.freqs),1,1]);  % target is: Mz = |alpha|^2-|beta|^2 = cos(theta)
        adj2.target_LFA=repmat(adj2.roi(adj2.ind) .* cos(tfa_rad),[length(opt.freqs),1]);
    else  % LS
        tmp = []; tmp2 = [];
        for i = 1:opt.nfreqs
            ttt = adj2.roi .* sin(tfa_rad) .* exp(1j*opt.tphase(:,:,i)); 
            tmp =[tmp; ttt];  % target is: Mt = 2*conj(alpha)*beta = sin(theta)*exp(1j*tphase)
            tmp2 = [tmp2; ttt(adj2.ind)]; 
        end
        adj2.targetImg_LFA = tmp;
        adj2.target_LFA = tmp2;
    end
else  % refocusing pulse
    if opt.do_mls==1  % MLS
        adj2.targetImg_LFA=repmat(double( adj2.roi>0 ),[length(opt.freqs), 1, 1]);  % target is: |beta|^2=1
        adj2.target_LFA = repmat(double( adj2.roi(adj2.ind)>0 ),[length(opt.freqs), 1]);
    else  % LS
        tmp = []; tmp2 = [];
        for i = 1:opt.nfreqs
            ttt=adj2.roi .* (-1.0) .* exp(1j*2*opt.tphase(:,:,i));  % target is: beta^2=-exp(1j*tphase)
            tmp = [tmp; ttt];
            tmp2 = [tmp2; ttt(adj2.ind)];
        end
        adj2.targetImg_LFA = tmp;
        adj2.target_LFA = tmp2;
    end
    %adj2.target_LFA=adj2.targetImg_LFA( adj2.ind );
end
%adj2.target_LFA=adj2.targetImg_LFA( adj2.ind );


% display B1 maps
nfigsx=ceil(sqrt(system.ncoils));
if nfigsx*(nfigsx-1)>system.ncoils
    nfigsy=nfigsx-1;
else
    nfigsy=nfigsx;
end
% if strcmpi(opt.display,'on')
%     to_disp=conv_3Dstack_2_2Dimg( permute(adj2.b1maps,[2 3 1]),nfigsx,nfigsy ) ;
%     
%     f=figure;
%     imagesc(abs(to_disp)); axis image; axis off; colormap(hot); colorbar;
%     title('Amplitude of B1 maps [T/V]');
% 
%     f=figure;
%     imagesc(angle(to_disp)); axis image; axis off; colormap(hot); colorbar; caxis([-pi pi]);
%     title('Phase of B1 maps [T/V]');
% end

% display B0 map
if strcmpi(opt.display,'on')
    figure;
    imagesc(adj2.b0map); axis image; axis off; colormap(hot); colorbar;
    title('B0 map [Hz]');
end

adj2.b1maps_=permute( adj2.b1maps(:,adj2.ind),[2 1] );
adj2.b0map_=adj2.b0map(adj2.ind);











