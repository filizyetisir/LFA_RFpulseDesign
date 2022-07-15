

function spokes=prepare_pulse(spokes,opt,system)
% function spokes=prepare_pulse(spokes,opt,system)

gamma=42.576*1e6;

spokes=prepare_gradient_spokes(spokes,system,opt);

spokes.sinc=spokes.sinc / max(abs(spokes.sinc));  % normalization of the RF amplitude
spokes.sumsinc=sum(spokes.sinc);
spokes.sumsincsq=sum(abs(spokes.sinc).^2);
spokes.sumsincabs=sum(abs(spokes.sinc));

if strcmpi(opt.display,'on')
    
    figure;
    times=[1:spokes.ntimes]'*system.deltat;
    subplot(2,1,1); plot(times*1e3,spokes.g(:,1)*1e3,'Color','blue','LineWidth',2); hold on;
    plot(times*1e3,spokes.g(:,2)*1e3,'Color','red','LineWidth',2);
    plot(times*1e3,spokes.g(:,3)*1e3,'Color','green','LineWidth',2);
    legend('Gx','Gy','Gz','Location','EastOutside'); legend('boxoff');
    title('Gradient trajectory [mT/m]');

    subplot(2,1,2); plot(times*1e3,spokes.k(:,1),'Color','blue','LineWidth',2); hold on;
    plot(times*1e3,spokes.k(:,2),'Color','red','LineWidth',2);
    plot(times*1e3,spokes.k(:,3),'Color','green','LineWidth',2);
    legend('kx','ky','kz','Location','EastOutside'); legend('boxoff');
    title('K-space trajectory [1/m]'); xlabel('Time [ms]');
    
%     figure;
%     plot(times(1:spokes.ntimesperspoke)*1e3,spokes.sinc,'LineWidth',2);
%     xlabel('Time [ms]');
%     title('Spoke sub-pulse profile');
%     
    zs=linspace(-2*spokes.sthick,2*spokes.sthick,500)';
    
    [mxf myf mzf]=bloch_simulation_z_correct(spokes,system,zs,0);
    MT1=abs( mxf + 1j*myf );    
    
    [mxf myf mzf]=bloch_simulation_z_correct(spokes,system,zs,50);
    MT2=abs( mxf + 1j*myf );

    [mxf myf mzf]=bloch_simulation_z_correct(spokes,system,zs,100);
    MT3=abs( mxf + 1j*myf );

    norm_fac=max(MT1);
    
    figure;
    plot(zs*1e3,MT1/norm_fac,'-r','LineWidth',2); hold on;
    plot(zs*1e3,MT2/norm_fac,'-g','LineWidth',2);
    plot(zs*1e3,MT3/norm_fac,'-b','LineWidth',2);
    
    xlabel('Position in z [mm]');
    title('Slice profile (STA regime)');
    
    box off;
    legend('0 Hz','50 Hz','100 Hz');
    
end


























