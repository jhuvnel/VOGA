MVI_path = '\\win.ad.jhu.edu\cloud\vnelhuman$\MVI\Study Subjects';
VOG_fdir = unique(extractfield(dir([MVI_path,filesep,'MVI*',filesep,'Visit*',filesep,'Rotary Chair',filesep,'*',filesep,'*LHRH-n240dps*']),'folder'));
for i = 1:length(VOG_fdir)
    Path = VOG_fdir{i};
    disp([num2str(i),'/',num2str(length(VOG_fdir)),': ',Path])    
    VOGA__RenameFiles('LHRH-n240dps','RH-240dps',Path);
    VOGA__RenameFiles('LHRH-240dps','LH-240dps',Path);
    VOGA__RenameFiles('LHRH-n120dps','RH-120dps',Path);
    VOGA__RenameFiles('LHRH-120dps','LH-120dps',Path);
    VOGA__RenameFiles('LHRH-n60dps','RH-60dps',Path);
    VOGA__RenameFiles('LHRH-60dps','LH-60dps',Path);
end
disp('Successfully finished!!')