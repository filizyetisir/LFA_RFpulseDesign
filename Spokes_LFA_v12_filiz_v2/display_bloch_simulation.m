

function display_bloch_simulation(MT,MZ,FA,BETA2,opt,adj,title_str)
% function display_bloch_simulation(MT,MZ,FA,BETA2,opt,adj,title_str)

if strcmpi('on',opt.display)==1

    
    % Whole ROI
%     
%     f=figure;
% 
%     subplot(3,2,1); 
%     imagesc(abs(MT)); axis image off; colormap hot; caxis([0 1]); colorbar; title(sprintf('|MT| (%s)',title_str));
% 
%     subplot(3,2,2); 
%     imagesc(angle(MT)); axis image off; colormap hot; caxis([-pi pi]); colorbar; title(sprintf('angle MT (%s)',title_str));
% 
%     subplot(3,2,3); 
%     imagesc(MZ); axis image off; colormap hot; caxis([-1 1]); colorbar; title(sprintf('|MZ| (%s)',title_str));
% 
%     subplot(3,2,4); 
%     if(opt.refoc_pulse == 0)
%         imagesc(FA); axis image off; colormap hot; caxis([0 opt.tfa*1.5]); colorbar; title(sprintf('FA (%s) [deg.]',title_str));
%     else
%         imagesc(FA); axis image off; colormap hot; caxis([0 180]); colorbar; title(sprintf('FA (%s) [deg.]',title_str));
%     end
%     subplot(3,2,5); 
%     imagesc(abs(BETA2)); axis image off; colormap hot; caxis([0 1]); colorbar; title(sprintf('|BETA^2| (%s)',title_str));
% 
%     subplot(3,2,6); 
%     imagesc(angle(BETA2)); axis image off; colormap hot; caxis([-pi pi]); colorbar; title(sprintf('angle BETA^2 (%s)',title_str));
%     
    
    
    % Designed ROI
    MT = MT.*adj.roi;
    MZ = MZ.*adj.roi;
    FA = FA.*adj.roi;
    BETA2 = BETA2.*adj.roi;
    
    f=figure;

    subplot(3,2,1); 
    imagesc(abs(MT)); axis image off; colormap hot; caxis([0 1]); colorbar; title(sprintf('|MT| (%s)',title_str));

    subplot(3,2,2); 
    imagesc(angle(MT)); axis image off; colormap hot; caxis([-pi pi]); colorbar; title(sprintf('angle MT (%s)',title_str));

    subplot(3,2,3); 
    imagesc(MZ); axis image off; colormap hot; caxis([-1 1]); colorbar; title(sprintf('|MZ| (%s)',title_str));

    subplot(3,2,4); 
    if(opt.refoc_pulse == 0)
        imagesc(FA); axis image off; colormap hot; caxis([0 opt.tfa*1.5]); colorbar; title(sprintf('FA (%s) [deg.]',title_str));
    else
        imagesc(FA); axis image off; colormap hot; caxis([0 180]); colorbar; title(sprintf('FA (%s) [deg.]',title_str));
    end
    subplot(3,2,5); 
    imagesc(abs(BETA2)); axis image off; colormap hot; caxis([0 1]); colorbar; title(sprintf('|BETA^2| (%s)',title_str));

    subplot(3,2,6); 
    imagesc(angle(BETA2)); axis image off; colormap hot; caxis([-pi pi]); colorbar; title(sprintf('angle BETA^2 (%s)',title_str));


end



