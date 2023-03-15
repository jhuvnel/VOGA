MVI_path = '\\win.ad.jhu.edu\cloud\vnelhuman$\MVI\Study Subjects';
VOG_fdir = unique(extractfield(dir([MVI_path,filesep,'MVI*',filesep,'Visit*',filesep,'eeVOR',filesep,'*',filesep,'*LHRH-n240dps*']),'folder'));
for i = 1:length(VOG_fdir)
    Path = VOG_fdir{i};
    disp([num2str(i),'/',num2str(length(VOG_fdir)),': ',Path])    
    VOGA__RenameFiles('LHRH-n240dps','RH-240dps',Path);
    VOGA__RenameFiles('LHRH-240dps','LH-240dps',Path);
end
disp('Successfully finished!!')