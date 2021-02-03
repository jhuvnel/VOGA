function [ang_mean, ang_std, ang_all] = calc_misalignment(stim_vect,V)
    if any(size(stim_vect)==1)
        M = repmat(stim_vect,size(V,1),1); %Don't invert
    else
        M = stim_vect; %Compare two vectors, do not invert
    end
    mag_V = sqrt((V(:,1).^2 + V(:,2).^2 + V(:,3).^2));
    mag_M = sqrt((M(:,1).^2 + M(:,2).^2 + M(:,3).^2));
    ang_all = acosd(dot(V,M,2)./(mag_V.*mag_M));
    ang_mean = mean(ang_all);
    ang_std = std(ang_all);
end