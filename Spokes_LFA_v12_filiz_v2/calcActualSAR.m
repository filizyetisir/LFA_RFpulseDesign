function sar = calcActualSAR()
%function sar = calcActualSAR()

% This function calculates the actual SAR from the SAR matrices not VOPs




%% INITIALIZE

nx = 100;
ny = 100;
nz = 100;

lx = 0.4;
ly = 0.4;
lz = 0.4;

nGram = 10; % 10 g average SAR

C = 8;
rows = 1;



%% POINT SAR

% Number of rows, coils, number of voxels in x,y,z, E field locations, ?
!/autofs/cluster/wald/6/filiz/Bay5_Utilities/ActualSARCalc/ComputeSAR2/bin/ComputeSAR2 1 8 pTXArbitrary_LTA 100 100 100 /autofs/cluster/wald/6/filiz/SimulatedDataSets/DataSet_100_100_100_7T_8Chn_Head/EField /autofs/cluster/wald/6/filiz/SimulatedDataSets/DataSet_100_100_100_7T_8Chn_Head/EField


sprintf('point sar map read successfully')


%% AVERAGE SAR

% Number of voxels in x y z, Length in x y z, N gram Average, density map
% file, sar map file

!/autofs/cluster/wald/6/filiz/Bay5_Utilities/ActualSARCalc/ComputeAverageSARMap2/bin/ComputeAverageSARMap 100 100 100 0.4 0.4 0.4 10 /autofs/cluster/wald/6/filiz/SimulatedDataSets/DataSet_100_100_100_7T_8Chn_Head/EField/density_100_100_100 sar_100_100_100 ./






%% GLOBAL SAR AND MAX LOCAL SAR AND SAR MAPS


pointSARMap = reshape(read_fdata('sar_100_100_100',-1), 100, 100, 100);
avgSARMap = reshape(read_fdata('avsar_100_100_100',-1), 100, 100, 100);

N = length(nonzeros(pointSARMap));
N2 = length(nonzeros(avgSARMap));

sar.maxLocSARav = max(max(max(avgSARMap)));
sar.maxLocSAR = max(max(max(pointSARMap)));

sar.globSARav = sum(sum(sum(avgSARMap)))/N2;
sar.globSAR = sum(sum(sum(pointSARMap)))/N;

sar.pointSARMap = pointSARMap;
sar.avgSARMap = avgSARMap;




end

