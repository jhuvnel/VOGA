%Make the general copy of this table
%Default table, no presets
trace_names = {'LX','RX','LY','RY','LZ','RZ','LLARP','RLARP','LRALP','RRALP','ALL'};
pos_filts = {'median','spline','sgolay1','sgolay2'}; %add more as needed
vel_filts = {'median','spline','sgolay1','sgolay2','irlssmooth'}; %add more as needed
filt1.pos = array2table(NaN(length(trace_names),length(pos_filts)));
filt1.pos.Properties.VariableNames = pos_filts;
filt1.pos.Properties.RowNames = trace_names;
filt1.vel = array2table(NaN(length(trace_names),length(vel_filts)));
filt1.vel.Properties.VariableNames = vel_filts;
filt1.vel.Properties.RowNames = trace_names;
%LDVOG best guess
LDVOG_filt1 = filt1;
LDVOG_filt1.pos.median = [11;11;3;3;3;3;NaN;NaN;NaN;NaN;NaN];
LDVOG_filt1.pos.spline = [1;1;0.999995;0.999995;0.9999995;0.9999995;NaN;NaN;NaN;NaN;NaN];
%NKI best guess
NKI_filt1 = filt1;
NKI_filt1.pos.median = [11;11;3;3;3;3;NaN;NaN;NaN;NaN;NaN];
NKI_filt1.pos.spline = [1;1;0.999995;0.999995;0.9999995;0.9999995;NaN;NaN;NaN;NaN;NaN];
%GNO best guess
GNO_filt1 = filt1;
GNO_filt1.vel.sgolay1(end) = 2;
GNO_filt1.vel.sgolay2(end) = 7;
GNO_filt1.vel.irlssmooth(end) = 5;
%Make filt_params struct
filt_params.default.filt1 = filt1;
filt_params.default.YLim.Pos = [];
filt_params.default.YLim.Vel = [];
filt_params.LDVOG.filt1 = LDVOG_filt1;
filt_params.LDVOG.YLim.Pos = [-30 30];
filt_params.LDVOG.YLim.Vel = [-50 50];
filt_params.NKI.filt1 = NKI_filt1;
filt_params.NKI.YLim.Pos = [-30 30];
filt_params.NKI.YLim.Vel = [-50 50];
filt_params.GNO.filt1 = GNO_filt1;
filt_params.GNO.YLim.Pos = [-30 30];
filt_params.GNO.YLim.Vel = [-250 250];
save([userpath,filesep,'VOGA_DefaultFilterParamsLocal.mat'],'filt_params')