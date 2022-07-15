

function A=compute_system_matrix(adj,system,spokes,freq_offset,SP_mean)
% function A=compute_system_matrix(adj,system,spokes,freq_offset,SP_mean)

gamma=42.576*1e6;
nchannels=system.ncoils;

A=zeros(adj.nnonzeropixels,spokes.nspokes*nchannels);
for chn=1:nchannels
    for sp=1:spokes.nspokes
        deltab0=( spokes.ntimes-spokes.sinc_time_and_spoke_to_time(:,sp) )*system.deltat*( adj.b0map(adj.ind')+freq_offset );
        ktrans=spokes.k( spokes.sinc_time_and_spoke_to_time(1,sp),1:2 );
        A(:,chn+(sp-1)*nchannels)=1j*2*pi*gamma*adj.b1maps(chn,adj.ind)*system.deltat ...
            .* exp( 1j*2*pi*(ktrans*adj.nonzeropixels(:,4:5)') ) ...
            .* ( spokes.sinc.' * exp( 1j*2*pi*deltab0 ) ) ...
            * SP_mean;
    end
end

