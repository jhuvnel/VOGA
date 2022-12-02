%% All ESC Notes
% Make Notes files/Segments for all EyeSeeCam folders
%% Make Notes Files/Segment
% ESC_dir = unique(extractfield(dir('\\win.ad.jhu.edu\cloud\vnelhuman$\MVI\Study Subjects\MVI*R*\Visit*\vHIT\ESC'),'folder'));
% for f = 1:length(ESC_dir)
%     Path = ESC_dir{f};
%     disp([num2str(f),'/',num2str(length(ESC_dir)),': ',Path])
%     %VOGA__makeFolders(Path,1,0); 
%     %flag = VOGA__ProcessRawData(Path);
%     %flag = MakeNotes([Path,filesep,'Raw Files']);
%     VOGA__Segment(Path,1);
% end
%% Find folders with Segments in them
ESC_dir = unique(extractfield(dir('\\win.ad.jhu.edu\cloud\vnelhuman$\MVI\Study Subjects\MVI*R*\Visit*\vHIT\ESC'),'folder'));
ESC_seg_dir = strrep(unique(extractfield(dir('\\win.ad.jhu.edu\cloud\vnelhuman$\MVI\Study Subjects\MVI*R*\Visit*\vHIT\ESC\Segmented Files\*.mat'),'folder')),'\Segmented Files','');
missing_segments = ESC_dir(~ismember(ESC_dir,ESC_seg_dir));
%% Rename .mat files with the .pdf names
% mat_files = extractfield(dir('*.mat'),'name');
% pdf_files = extractfield(dir('*.pdf'),'name');
% date_pat1 = digitsPattern(4)+'-'+digitsPattern(2)+'-'+digitsPattern(2)+'_'+digitsPattern(2)+'.'+digitsPattern(2)+'.'+digitsPattern(2);
% for j = 1:length(mat_files)
%     fname = mat_files{j};
%     date = extract(fname,date_pat1);
%     rel_pdf = pdf_files(contains(pdf_files,date));
%     if ~isempty(rel_pdf)
%         movefile(fname,strrep(rel_pdf{:},'.pdf','_export.mat'))
%     end
% end
%% Replace SettingA/B with Setting
% VOGA__RenameFiles('Setting A','MotionMod')
% VOGA__RenameFiles('Setting B','ConstantRate')
%% Add xls extension
% file_Path = extractfield(dir,'name',~extractfield(dir,'isdir'));
% file_Path(contains(file_Path,{'export','.pdf','DS_Store','.txt'})) = [];
% for j = 1:length(file_Path)
%     rel_file = file_Path{j};
%     movefile([cd,filesep,rel_file],[cd,filesep,rel_file,'.xls'])
% end

% ADD mat extension
% file_Path = extractfield(dir,'name',~extractfield(dir,'isdir'));
% file_Path(contains(file_Path,{'.xls','.pdf','DS_Store','.mat','.txt'})) = [];
% for j = 1:length(file_Path)
%     rel_file = file_Path{j};
%     movefile([cd,filesep,rel_file],[cd,filesep,rel_file,'.mat'])
% end

%Remove .mat extension
% file_Path = extractfield(dir('*.mat'),'name');
% for j = 1:length(file_Path)
%     rel_file = file_Path{j};
%     movefile([cd,filesep,rel_file],[cd,filesep,strrep(rel_file,'.mat','')])
% end
%% Rename Folders Titled EyeSeeCam as ESC
% vHIT_dir = dir('\\win.ad.jhu.edu\cloud\vnelhuman$\MVI\Study Subjects\MVI*R*\Visit*\vHIT');
% vHIT_dir(contains(extractfield(vHIT_dir,'name'),{'GNO','.'})) = [];
% vHIT_file = strcat(extractfield(vHIT_dir,'folder'),filesep,extractfield(vHIT_dir,'name'));
% old_name = vHIT_file(contains(vHIT_file,'EyeSeeCam'));
% for i = 1:length(old_name)
%     movefile(old_name{i},strrep(old_name{i},'EyeSeeCam','ESC'))
% end