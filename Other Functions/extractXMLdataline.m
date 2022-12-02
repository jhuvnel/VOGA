function dat = extractXMLdataline(line)
    if iscell(line)
        line = line{:};
    end
    l_c = strfind(line,'<');
    r_c = strfind(line,'>');
    dat = line((r_c(1)+1):(l_c(2)-1));    
end