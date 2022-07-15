


function stop=display_iteration(x,optimValues,state)
% function stop=display_iteration(x,optimValues,state)

stop=false;

nchannels=size(opt.b1maps,1);
rf=x(1:nchannels*spokes.nspokes,1);  % spokes coeficients are stored first
spokes.spcoords=reshape( x(nchannels*spokes.nspokes+1:end),spokes.nspokes,3 );  % spokes coordinates are stored after

% compute gradient trajectory corresponding to these spokes coordinates
prepare_gradient(spokes,system);

% compute system matrix
A=compute_system_matrix(opt,system,spokes);

% compute achieved magnetization
mt=A*rf;
mtimg=zeros(opt.nx,opt.ny,opt.nz);
for n=1:opt.nnonzeropixels
    [i j k x y z]=get_coord_non_zero_pixel(opt,n);
    mtimg(i,j,k)=mt(n);
end

figure; cmax=max(max(abs(opt.mttimg)))*1.5;
subplot(2,2,1); imagesc(abs(opt.mttimg)); colormap(hot); colorbar; axis image; caxis([0 cmax]); axis off; title('Target magnetization (mag.)');
subplot(2,2,3); imagesc(angle(opt.mttimg)); colormap(hot); colorbar; axis image; caxis([-pi pi]); axis off; title('Target magnetization (pha.)');
subplot(2,2,2); imagesc(abs(mtimg)); colormap(hot); colorbar; axis image; caxis([0 cmax]); axis off; title('Achieved magnetization (mag.)');
subplot(2,2,4); imagesc(angle(mtimg)); colormap(hot); colorbar; axis image; caxis([-pi pi]); axis off; title('Achieved magnetization (pha.)');







