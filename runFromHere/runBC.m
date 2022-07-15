
center = 0;
avg = 0;
avg_noverse = 1;

verse_factor = 1;
display = 1;

%% BC center
if(center == 1)
    % 180
    spokesdefpath = 'spokes_def_180.txt';
    BCcenter_180(verse_factor, display, spokesdefpath)
    mkdir ./Results/BCcenter_BCcenter/180
    load all.mat
    save('./Results/BCcenter_BCcenter/180/all', 'pulse', 'opt', 'adj', 'system', 'spokes', 'errorr');
    
    % 90
    spokesdefpath = 'spokes_def_90.txt';
    BCcenter_90(verse_factor, display, spokesdefpath)
    mkdir ./Results/BCcenter_BCcenter/90
    load all.mat
    save('./Results/BCcenter_BCcenter/90/all', 'pulse', 'opt', 'adj', 'system', 'spokes', 'errorr');
end

%% BC avg
if(avg == 1)
    % 180
    spokesdefpath = 'spokes_def_180.txt';
    BCavg_180(verse_factor, display, spokesdefpath)
    mkdir ./Results/BCavg_BCavg/180
    load all.mat
    save('./Results/BCavg_BCavg/180/all', 'pulse', 'opt', 'adj', 'system', 'spokes', 'errorr');
    
    % 90
    spokesdefpath = 'spokes_def_90.txt';
    BCavg_90(verse_factor, display, spokesdefpath)
    mkdir ./Results/BCavg_BCavg/90
    load all.mat
    save('./Results/BCavg_BCavg/90/all', 'pulse', 'opt', 'adj', 'system', 'spokes', 'errorr');
end

%% BC avg no VERSE
if(avg_noverse == 1)
    verse_factor = -1;
    % 180
    spokesdefpath = 'spokes_def_180_noverse.txt';
    BCavg_180(verse_factor, display, spokesdefpath)
    mkdir ./Results/noverse_BCavg_BCavg/180
    load all.mat
    save('./Results/noverse_BCavg_BCavg/180/all', 'pulse', 'opt', 'adj', 'system', 'spokes', 'errorr');
    
    % 90
    spokesdefpath = 'spokes_def_90_noverse.txt';
    BCavg_90(verse_factor, display, spokesdefpath)
    mkdir ./Results/noverse_BCavg_BCavg/90
    load all.mat
    save('./Results/noverse_BCavg_BCavg/90/all', 'pulse', 'opt', 'adj', 'system', 'spokes', 'errorr');
end