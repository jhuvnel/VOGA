%% Make GNO Summary Table
% Based on GNO_XML_Creator_v3 and GNO_SummaryTableMaker_v5 in the vHIT analysis scripts
function [tab,Data] = MakeGNOSummaryTable(out_path,XML_Path,rerun)
% Handle inputs
if nargin < 3
  rerun = 1; 
end
if nargin < 1
    out_path = cd;
end
if nargin < 2 || isempty(XML_Path)
    XML_Path = [out_path,filesep,'GNO',filesep,'Raw Files'];
end
if isfile([out_path,filesep,'GNO_SystemOutput.mat'])&&isfile([out_path,filesep,'GNO_Summary.mat'])&&~rerun
    load([out_path,filesep,'GNO_SystemOutput.mat'],'Data')
    load([out_path,filesep,'GNO_Summary.mat'],'tab')
else
    [tab,Data] = GNO_XML_Parser(XML_Path);
    if ~isempty(tab)
        save([out_path,filesep,'GNO_SystemOutput.mat'],'Data')
        save([out_path,filesep,'GNO_Summary.mat'],'tab')
    end
end
end