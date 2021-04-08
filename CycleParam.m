%% CycleParam.m 
%This function takes in a CycAvg file called by MakeCycAvg.m
%and parameterizes each cycle, adding results to CycAvg struct.

%Impulse Fit (Many params including gain, Latency, Saccade Analysis)
%Sine Fit (Gain, Phase, Misalignment)
%Exponential Fit (Tau, Magnititude)
%Magnitude (Max Magnitude, Misalignment)
function CycAvg = CycleParam(CycAvg)
    %Figure out what kind of file it is and what should be done
    if contains(CycAvg.name,'Impulse') %Impulses with many params
        if contains(CycAvg.info.goggle_ver,'GNO') %Only 1D eye tracking really
            
        else %LDVOG, or NKI are in 3D
            
        end
    elseif contains(CycAvg.name,'Sin') %Sine Fits
        
    elseif contains(CycAvg.name,{'Step','Activation'}) %Velocity Step/Activation
        
    else %Pulse Trains, 65 Vector
        
    end
end
