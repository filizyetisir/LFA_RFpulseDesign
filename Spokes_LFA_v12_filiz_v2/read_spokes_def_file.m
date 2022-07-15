
function [adj sys optimization spokes sar]=read_spokes_def_file(filepath, optimization)
% function [adj sys optimization spokes sar]=read_spokes_def_file(filepath)

file=fopen(filepath,'r');

%% Sys parameters

tline=fgets(file);
tline=fgets(file); sys.gmax=str2num( tline(26:size(tline,2)-2) );
sys.gmax=sys.gmax*1e-3;  % mT/m -> T/m

tline=fgets(file); sys.smax=str2num( tline(28:size(tline,2)-2) );
tline=fgets(file); sys.umax=str2num( tline(38:size(tline,2)-2) );
sys.umax=sys.umax;
tline=fgets(file); sys.ncoils=str2num( tline(20:size(tline,2)-2) );
tline=fgets(file); sys.deltat=tline_to_param(tline,19,size(tline,2)-1);
tline=fgets(file); sys.delay_rf_grad=tline_to_param(tline,36,size(tline,2)-1);



%% Adjustment data

tline=fgets(file); tline=fgets(file);

%whole data
tline=fgets(file); wholepath=tline(17:size(tline,2)-2) ;
eval(sprintf('load %s',wholepath));
downsample = optimization.downsample;
dsr = optimization.dsr; % downsample rate

if(downsample)
    adj.roiWhole = roiWhole(1:dsr:end,1:dsr:end);
    adj.b0mapWhole = b0mapWhole(1:dsr:end,1:dsr:end);
    tmp = b1mapsWhole(1:dsr:end,1:dsr:end,:,:);
    adj.b1mapsWhole=( permute( tmp,[4 1 2 3] ) );
else
    adj.roiWhole = roiWhole;
    adj.b0mapWhole = b0mapWhole;
    adj.b1mapsWhole=( permute( b1mapsWhole,[4 1 2 3] ) );
end

%inner brain data
tline=fgets(file); innerpath=tline(28:size(tline,2)-2) ;
eval(sprintf('load %s',innerpath));

if(downsample)
    adj.roiInner = roiInner(1:dsr:end,1:dsr:end);
    adj.b0mapInner = b0mapInner(1:dsr:end,1:dsr:end);
    tmp = b1mapsInner(1:dsr:end,1:dsr:end,:,:);
    adj.b1mapsInner=( permute( tmp,[4 1 2 3] ) );
else
    adj.roiInner = roiInner;
    adj.b0mapInner = b0mapInner;
    adj.b1mapsInner=( permute( b1mapsInner,[4 1 2 3] ) );
end

%spcial region inside brain data
tline=fgets(file); specialpath=tline(39:size(tline,2)-2) ;
eval(sprintf('load %s',specialpath));

if(downsample)
    adj.roiSpecial = roiSpecial(1:dsr:end,1:dsr:end);
    adj.b0mapSpecial = b0mapSpecial(1:dsr:end,1:dsr:end);
    tmp = b1mapsSpecial(1:dsr:end,1:dsr:end,:,:);
    adj.b1mapsSpecial=( permute( tmp,[4 1 2 3] ) );
else
    adj.roiSpecial = roiSpecial;
    adj.b0mapSpecial = b0mapSpecial;
    adj.b1mapsSpecial=( permute( b1mapsSpecial,[4 1 2 3] ) );
end

% Check
if size(adj.b1mapsWhole,1)~=sys.ncoils
    error('Inconsistent number of channels in B1+ maps.');
end

% x,y,z dimensions in voxel
adj.nx=size(adj.b1mapsWhole,2);
adj.ny=size(adj.b1mapsWhole,3);
adj.nz=1;

% SODA
adj.SODA=SODA;

% Position Info
if(downsample)
    adj.xs=SODA.rows_pos(1:dsr:end)*0.001;  % mm -> m. Also, in this code x is the row index, whereas in the SODA file x is x in the dicom convention.
    % Therfore x and y have to be swapped.
    adj.ys=SODA.cols_pos(1:dsr:end)*0.001;  % mm -> m
    adj.zs=0.0;  % always design pulse @ z=0
    
else
    adj.xs=SODA.rows_pos*0.001;  % mm -> m. Also, in this code x is the row index, whereas in the SODA file x is x in the dicom convention.
    % Therfore x and y have to be swapped.
    adj.ys=SODA.cols_pos*0.001;  % mm -> m
    adj.zs=0.0;  % always design pulse @ z=0
end

% Pixel width
adj.dx=abs( adj.xs(2)-adj.xs(1) );
adj.dy=abs( adj.ys(2)-adj.ys(1) );
adj.dz=0.005;  % place-holder

% x,y,z dimensions in m
if(downsample)
    adj.lx=abs( adj.xs(end)-adj.xs(1) ) + adj.dx/dsr;
    adj.ly=abs( adj.ys(end)-adj.ys(1) ) + adj.dy/dsr;
    adj.lz=abs( adj.zs(end)-adj.zs(1) ) + adj.dz/dsr;
else
    adj.lx=abs( adj.xs(end)-adj.xs(1) ) + adj.dx;
    adj.ly=abs( adj.ys(end)-adj.ys(1) ) + adj.dy;
    adj.lz=abs( adj.zs(end)-adj.zs(1) ) + adj.dz;
end


%% Pulse Parameters

tline=fgets(file); optimization.tfa=tline_to_param(tline,19,size(tline,2)-2);

% spokes parameters
tline=fgets(file); tline=fgets(file);
tline=fgets(file); spokes.nspokes=tline_to_param(tline,18,size(tline,2)-2);
spokes.sdir=zeros(3,1);
tline=fgets(file); [spokes.sdir(1) spokes.sdir(2) spokes.sdir(3)]=tline_to_3params(tline,17,size(tline,2)-2);
tline=fgets(file); spokes.sthick=tline_to_param(tline,20,size(tline,2)-2);
tline=fgets(file); spokes.soff=tline_to_param(tline,34,size(tline,2)-2);
tline=fgets(file); spokes.nsinczc=tline_to_param(tline,31,size(tline,2)-2);
tline=fgets(file); spokes.SS_grad_down=tline_to_param(tline,58,size(tline,2)-2)/100;
spokes.spcoords=zeros(spokes.nspokes,3);
for i=1:spokes.nspokes
    tline=fgets(file); [spokes.spcoords(i,1) spokes.spcoords(i,2) spokes.spcoords(i,3)]=tline_to_3params(tline,22,size(tline,2)-2);
end



%% Opt. parameters

tline=fgets(file); tline=fgets(file);
tline=fgets(file); optimization.niters=tline_to_param(tline,30,size(tline,2)-2);



%% SAR control

tline=fgets(file); tline=fgets(file);
tline=fgets(file); sar.lsarpath=tline(36:size(tline,2)-2);
tline=fgets(file); sar.lsarmax=tline_to_param(tline,36,size(tline,2)-2);
tline=fgets(file); sar.gsarpath=tline(35:size(tline,2)-2);
tline=fgets(file); sar.gsarmax=tline_to_param(tline,37,size(tline,2)-2);
tline=fgets(file); sar.ppmax=tline_to_param(tline,24,size(tline,2)-2);
tline=fgets(file); sar.ppavmax=tline_to_param(tline,32,size(tline,2)-2);


%% Results

tline=fgets(file); tline=fgets(file);
tline=fgets(file); optimization.resdir=tline(18:size(tline,2)-2);

fclose(file);










function [B1, B2, B3] = transform_coordinates(A1, A2, A3, Ma2b)
% function [B1, B2, B3] = transform_coordinates(A1, A2, A3, Ma2b)

B1 = Ma2b(1,1)*A1 + Ma2b(1,2)*A2 + Ma2b(1,3)*A3 + Ma2b(1,4);
B2 = Ma2b(2,1)*A1 + Ma2b(2,2)*A2 + Ma2b(2,3)*A3 + Ma2b(2,4);
B3 = Ma2b(3,1)*A1 + Ma2b(3,2)*A2 + Ma2b(3,3)*A3 + Ma2b(3,4);













