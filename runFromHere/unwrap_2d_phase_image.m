

function [img_unwrapped]=unwrap_2d_phase_image(img,mask,method,refpoint)
% function [img_unwrapped]=unwrap_2d_phase_image(img,mask,method,refpoint)
% This function unwrap a 2D phase image using either Golstein algorithm or
% the quality guided approach.
%
% Inputs:
% img: phase image to unwrap
% mask: mask within which to perform the phase unwrapping
% method: select 1 for Golstein's approach and any other value for the
% quality guided approach
% refpoint: select 1 for using the center of the phase image as the
% reference phase point
%
% Outputs:
% img_unwrapped: the unwrap phase image

max_box_radius=4;
threshold_std=5;

if method==1
    % Golstein approach
    residue_charge=PhaseResidues(img,mask);                            % calculate phase residues
    branch_cuts=BranchCuts(residue_charge,max_box_radius,mask);        % place branch cuts
    [img_unwrapped,rowref,colref]=FloodFill(img,branch_cuts,mask,refpoint);     % flood fill phase unwrapping
    
else
    % quality guided approach
    img_unwrapped=zeros(size(img));                                    % zero starting matrix for unwrapped phase
    adjoin=zeros(size(img));                                           % zero starting matrix for adjoin matrix
    unwrapped_binary=zeros(size(img));                                 % binary image to mark unwrapped pixels
    
    im_phase_quality=PhaseDerivativeVariance(img);   
    minp=im_phase_quality(2:end-1, 2:end-1); minp=min(minp(:));
    maxp=im_phase_quality(2:end-1, 2:end-1); maxp=max(maxp(:));
    if(refpoint)
        figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), axis square, axis off; title('Phase quality map'); 
        uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
        [xpoint,ypoint] = ginput(1);                                       % select starting point for the guided floodfill algorithm
    else
        xpoint=size(img,1)/2;
        ypoint=size(img,2)/2;
    end
        
    colref=round(xpoint); rowref=round(ypoint);
    img_unwrapped(rowref,colref)=img(rowref,colref);                   % save the unwrapped values
    unwrapped_binary(rowref,colref,1)=1;
    if mask(rowref-1, colref, 1)==1 adjoin(rowref-1, colref, 1)=1; end       % mark the pixels adjoining the selected point
    if mask(rowref+1, colref, 1)==1 adjoin(rowref+1, colref, 1)=1; end
    if mask(rowref, colref-1, 1)==1 adjoin(rowref, colref-1, 1)=1; end
    if mask(rowref, colref+1, 1)==1 adjoin(rowref, colref+1, 1)=1; end
    img_unwrapped=GuidedFloodFill(img,img_unwrapped,unwrapped_binary,im_phase_quality,adjoin,mask);    % unwrap
        
end
