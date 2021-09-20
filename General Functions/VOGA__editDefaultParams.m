function VOGA__editDefaultParams
if ~any(contains(extractfield(dir(userpath),'name'),'VOGA_DefaultFilterParamsLocal.mat'))
    generateFiltParams;
end
load('VOGA_DefaultFilterParamsLocal.mat','filt_params')
fields = fieldnames(filt_params);
fields(contains(fields,'default')) = []; %Don't edit the deafult
opts = [strcat({'Edit '},fields);{'Add New Goggle';'Save'}];
[ind,tf1] = nmlistdlg('PromptString','Select an action:',...
                   'SelectionMode','single',...
                   'ListSize',[150 125],...
                   'ListString',opts);  
while tf1
    if ind <= length(fields) %Edit goggle params
        opts2 = {'Filt Pos','Filt Vel','YLim'};
        [ind2,tf2] = nmlistdlg('PromptString','Select an action:',...
                   'SelectionMode','single',...
                   'ListSize',[150 125],...
                   'ListString',opts2); 
        if tf2
            if ind2 == 1 %Edit Position
                 prompt = {['Median',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Spline',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Sav-Gol 1',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
             ['Sav-Gol 2',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL'};
                dlgtitle = 'Filter position';
                definput = strrep(cellfun(@(x) num2str(x,10),table2cell(filt_params.(fields{ind}).filt1.pos),'UniformOutput',false),'NaN','');
                temp_filt_params_p = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 10],definput,'on',4,[0.5 0.5 3.5 7]));
                if ~isempty(temp_filt_params_p) %Didn't hit cancel
                    filt_params.(fields{ind}).filt1.pos{:,:} = reshape(temp_filt_params_p,11,[]);
                end
            elseif ind2 ==2 %Edit Vel
                prompt = {['Median',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Spline',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Sav-Gol 1',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
             ['Sav-Gol 2',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Irlssmooth',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL'};
                dlgtitle = 'Filter velocity';
                definput = strrep(cellfun(@(x) num2str(x,10),table2cell(filt_params.(fields{ind}).filt1.vel),'UniformOutput',false),'NaN','');
                temp_filt_params_v = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 10],definput,'on',5,[0.5 0.5 4.25 7]));
                if ~isempty(temp_filt_params_v) %Didn't hit cancel
                    filt_params.(fields{ind}).filt1.vel{:,:} = reshape(temp_filt_params_v,11,[]);
                end
            else %YLim
                prompt = {['Set Y-axis limits',newline,newline,'Position:',newline,newline,'Lower Limit:'],...
            'Upper Limit:',['Velocity:',newline,newline,'Lower Limit:'],'Upper Limit:'};
                dlgtitle = 'Y-axis Limits';
                definput = cellfun(@(x) num2str(x,10),num2cell([filt_params.(fields{ind}).YLim.Pos,filt_params.(fields{ind}).YLim.Vel]),'UniformOutput',false);
                out_nums = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 18],definput,'on',2,[0.5 0.5 3 2.25]));
                if ~isempty(out_nums)
                    filt_params.(fields{ind}).YLim.Pos = sort([out_nums(1),out_nums(2)]);
                    filt_params.(fields{ind}).YLim.Vel = sort([out_nums(3),out_nums(4)]);
                end
            end
        end
    elseif ind == length(fields)+1 %Add new 
        goggle_name = inputdlg('Name of the goggle set (short and capitalized)');
        if ~isempty(goggle_name)&&~isempty(goggle_name{:})
            gog = goggle_name{:};
            %Position Filter
                prompt = {['Median',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Spline',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Sav-Gol 1',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
             ['Sav-Gol 2',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL'};
                dlgtitle = 'Filter position';
                definput = strrep(cellfun(@(x) num2str(x,10),table2cell(filt_params.default.filt1.pos),'UniformOutput',false),'NaN','');
                temp_filt_params_p = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 10],definput,'on',4,[0.5 0.5 3.5 7]));
                if ~isempty(temp_filt_params_p) %Didn't hit cancel
                    filt_params.(gog).filt1.pos{:,:} = reshape(temp_filt_params_p,11,[]);
                else
                    filt_params.(gog).filt1.pos = file_params.deafult.filt1.pos;
                end
                %Velocity Filter
                prompt = {['Median',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Spline',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Sav-Gol 1',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
             ['Sav-Gol 2',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL',...
            ['Irlssmooth',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:','L L:','R L:','L R:','R R:','ALL'};
                dlgtitle = 'Filter velocity';
                definput = strrep(cellfun(@(x) num2str(x,10),table2cell(filt_params.default.filt1.vel),'UniformOutput',false),'NaN','');
                temp_filt_params_v = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 10],definput,'on',5,[0.5 0.5 4.25 7]));
                if ~isempty(temp_filt_params_v) %Didn't hit cancel
                    filt_params.(gog).filt1.vel{:,:} = reshape(temp_filt_params_v,11,[]);
                else
                    filt_params.(gog).filt1.vel = file_params.deafult.filt1.vel;
                end
                %Set YLim
                prompt = {['Set Y-axis limits',newline,newline,'Position:',newline,newline,'Lower Limit:'],...
            'Upper Limit:',['Velocity:',newline,newline,'Lower Limit:'],'Upper Limit:'};
                dlgtitle = 'Y-axis Limits';
                definput = cellfun(@(x) num2str(x,10),num2cell([-30,30,-100,100]),'UniformOutput',false);
                out_nums = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 18],definput,'on',2,[0.5 0.5 3 2.25]));
                if ~isempty(out_nums)
                    filt_params.(gog).YLim.Pos = sort([out_nums(1),out_nums(2)]);
                    filt_params.(gog).YLim.Vel = sort([out_nums(3),out_nums(4)]);
                end
        end
    elseif ind == length(opts)
        save([userpath,filesep,'VOGA_DefaultFilterParamsLocal.mat'],'filt_params')
        disp('VOGA_DefaultFilterParamsLocal.mat was saved')
    end
    %% Poll for new reponse
    fields = fieldnames(filt_params);
    fields(contains(fields,'default')) = []; %Don't edit the deafult
    opts = [strcat({'Edit '},fields);{'Add New Goggle';'Save'}];
    [ind,tf1] = nmlistdlg('PromptString','Select an action:',...
                   'SelectionMode','single',...
                   'ListSize',[150 125],...
                   'ListString',opts);    
end    
end