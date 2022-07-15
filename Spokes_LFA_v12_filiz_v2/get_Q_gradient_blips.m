

function a_gblips=get_Q_gradient_blips(spokes,adj,system)
% function a_gblips=get_Q_gradient_blips(spokes,adj,system)
% get the Q matrices (=alpha parameters since beta=0 for such blips)
% associated with the gradient blips used to move the spoke location. These
% are needed for computation of the Hessian. The method below is
% an explicit calculation and accounts for the impact of B0 on actual spoke
% positions

gamma=42.576*1e6;

% make sure this is complex, otherwise this creates error in the MEX function...
a_gblips=zeros(adj.nnonzeropixels,spokes.nspokes);

for i=1:spokes.nspokes
    
    % determine start and end time of G blip for current spoke
    timestart=spokes.sinc_time_and_spoke_to_time(end,i)+1;
    if i==spokes.nspokes
        timeend=spokes.ntimes;
    else
        timeend=spokes.sinc_time_and_spoke_to_time(1,i+1)-1;
    end
    
    % compute alpha parameter for this blip
    phi=0;
    for j=timestart:timeend        
        grad=adj.nonzeropixels(:,4:6) * spokes.g(j,:).' - adj.b0map_/gamma;
        phi=phi - 2*pi*gamma*system.deltat*grad;    
    end
    
    a_gblips(:,i) = exp( -1j*phi/2 );
end

a_gblips = complex( real(a_gblips),imag(a_gblips) );






