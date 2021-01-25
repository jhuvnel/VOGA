function keep_tr = MakeCycAvg__selectCycles(type,keep_tr)
%function [keep_tr,ha] = MakeCycAvg__selectCycles(ha,type,colors,line_wid,te,t_snip,stims,Fs,Data,info,filt,LE_V,RE_V,keep_inds,keep_tr,In_FileName)
    if type==1   
        cyc_num = 1:length(keep_tr);
        [ind,tf] = nmlistdlg('PromptString','Select cycles:',...
                           'SelectionMode','multiple',...
                           'InitialValue',cyc_num(keep_tr),...
                           'ListSize',[100 250],...
                           'ListString',cellstr(num2str(cyc_num')),...
                           'Position',[11,5.25,2,5.25]); 
        if tf==0
            return;
        else
            keep_tr = false(1,length(keep_tr));
            keep_tr(ind) = true;
        end       
        %% OLD VERSION WITH BUTTON CYCLE SELECTION
%         selection_type = questdlg('Cycle selection with buttons or keyboard shortcuts?','','Buttons','Keyboard','Buttons');
%         actions = {'Next Cycle','Previous Cycle','Remove Cycle','Add Cycle','Done'};
%         %Start with 1st cycle
%         cyc = 1;
%         if keep_tr(cyc)
%             set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_k)
%         else
%             set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_r)
%         end
%         set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.bold,'visible','on')
%         %Prompt
%         if strcmp(selection_type,'Buttons')
%             %Start with selecting things on and off on the graph
%             [ind,tf] = nmlistdlg('PromptString','Select an action:',...
%                            'SelectionMode','single',...
%                            'ListSize',[100 100],...
%                            'ListString',actions,...
%                            'Position',[11,7.75,2,2.75]); 
%             if tf==0
%                 command = actions{end};
%             else
%                 command = actions{ind};
%             end           
%         else %Keyboard shortcuts
%             clc; %clear any previous attempts
%             disp('Allowed commands:')
%             disp('n = next cycle, no change to cycle status')
%             disp('v = previous cycle, no change to cycle status')
%             disp('j = remove cycle')
%             disp('g = add cycle')
%             disp('d = done')
%             command = input('','s');
%         end
%         while ~(strcmp(command,'d')||strcmp(command,'Done'))
%             switch command
%                 case {'n','Next Cycle'}
%                     %Remove previous bold
%                     if keep_tr(cyc)
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_keep)
%                         set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.norm)
%                     else
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_rm)
%                         set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.norm)
%                         set(findobj(ha(4),'type','line','Tag',['Cycle_',num2str(cyc)]),'visible','off')
%                     end
%                     %Move cycle
%                     if cyc ~= length(keep_tr)
%                         cyc = cyc + 1;
%                     else
%                         cyc = 1; %Go back around the beginning
%                     end
%                     %Bold next cycle
%                     if keep_tr(cyc)
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_k)
%                     else
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_r)
%                     end
%                     set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.bold,'visible','on')
%                 case {'v','Previous Cycle'}
%                     %Remove previous bold
%                     if keep_tr(cyc)
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_keep)
%                         set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.norm)
%                     else
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_rm)
%                         set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.norm)
%                         set(findobj(ha(4),'type','line','Tag',['Cycle_',num2str(cyc)]),'visible','off')
%                     end
%                     if cyc ~= 1
%                         cyc = cyc - 1;
%                     else
%                         cyc = length(keep_tr); %Go to the end
%                     end
%                     %Bold next cycle
%                     if keep_tr(cyc)
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_k)
%                     else
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_r)
%                     end
%                     set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.bold,'visible','on')
%                 case {'j','Remove Cycle'}
%                     keep_tr(cyc) = false;
%                     %Remove bold
%                     set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_rm)
%                     set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.norm)
%                     set(findobj(ha(4),'type','line','Tag',['Cycle_',num2str(cyc)]),'visible','off')
%                     %Go to next cycle
%                     if cyc ~= length(keep_tr)
%                         cyc = cyc + 1;
%                     else
%                         cyc = 1; %Go back around the beginning
%                     end
%                     %Bold next cycle
%                     if keep_tr(cyc)
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_k)
%                     else
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_r)
%                     end
%                     set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.bold,'visible','on')
%                     %Recalc CycAvg
%                     title(ha(4),['Accepted Cycles: ',num2str(sum(keep_tr)),' of ',num2str(length(keep_tr))])
%                     CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
%                     MakeCycAvg__plotCycAvg(ha(5),type,colors,CycAvg);
%                 case {'g','Add Cycle'}
%                     keep_tr(cyc) = true;
%                     %Remove bold
%                     set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_keep)
%                     set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.norm)
%                     set(findobj(ha(4),'type','line','Tag',['Cycle_',num2str(cyc)]),'visible','on')
%                     %Go to next cycle
%                     if cyc ~= length(keep_tr)
%                         cyc = cyc + 1;
%                     else
%                         cyc = 1; %Go back around the beginning
%                     end
%                     %Bold next cycle
%                     if keep_tr(cyc)
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_k)
%                     else
%                         set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_bold_r)
%                     end
%                     set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.bold,'visible','on')
%                     %Recalc CycAvg and change title
%                     title(ha(4),['Accepted Cycles: ',num2str(sum(keep_tr)),' of ',num2str(length(keep_tr))])
%                     CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
%                     MakeCycAvg__plotCycAvg(ha(5),type,colors,CycAvg);
%                 otherwise
%                     disp('Input not recognized')
%             end
%             %Prompt
%             if strcmp(selection_type,'Buttons')
%                 %Start with selecting things on and off on the graph
%                 [ind,tf] = nmlistdlg('PromptString','Select an action:',...
%                                'SelectionMode','single',...
%                                'ListSize',[100 100],...
%                                'ListString',actions,...
%                                'Position',[11,7.75,2,2.75]); 
%                 if tf==0
%                     command = actions{end};
%                 else
%                     command = actions{ind};
%                 end           
%             else %Keyboard shortcuts
%                 command = input('','s');
%             end
%         end     
%         %Remove any bolded traces
%         if keep_tr(cyc)
%             set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_keep)
%             set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.norm)
%         else
%             set(findobj(gcf,'type','patch','Tag',['Cycle_',num2str(cyc)]),'FaceColor',colors.cyc_rm)
%             set(findobj(gcf,'type','line','Tag',['Cycle_',num2str(cyc)]),'LineWidth',line_wid.norm)
%             set(findobj(ha(4),'type','line','Tag',['Cycle_',num2str(cyc)]),'visible','off')
%         end
    else
        disp('No cycle selection available for this data type')
        keep_tr = true;
    end
end