function p = sinefitbounds(params)
%Takes in a 3x1 vector with half cycle fits and phase.
% Makes sure all gain amplitude values are positive and that the phase is
% adjusted accordingly
    if params(1) < 0 && params(2) < 0
        params(1) = abs(params(1));
        params(2) = abs(params(2));
        if params(3) > 0
            params(3) = params(3) - 180;
        else
            params(3) = params(3) + 180;
        end
    elseif params(1) < 0
        if abs(params(1)) > abs(params(2))
            params(1) = abs(params(1)) - params(2);
            params(2) = 0;
            if params(3) > 0
                params(3) = params(3) - 180;
            else
                params(3) = params(3) + 180;
            end
        else
            params(2) = params(2) - abs(params(1));
            params(1) = 0;
        end
    elseif params(2) < 0
        if abs(params(1)) > abs(params(2))
            params(1) = params(1) - abs(params(2));
            params(2) = 0;
        else
            params(2) = abs(params(2)) - params(1);
            params(1) = 0;
            if params(3) > 0
                params(3) = params(3) - 180;
            else
                params(3) = params(3) + 180;
            end
        end
    end
    %Now fix phase to be within -180 to 180
    while params(3) > 180
        params(3) = params(3) - 360;
    end
    while params(3) < -180
        params(3) = params(3) + 360;
    end
    p = params;
end

