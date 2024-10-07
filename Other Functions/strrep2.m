function str_out = strrep2(str_in,strs,str)
if nargin < 3
   error('Not enough inputs.')
end
str_out = str_in;
if iscell(strs)
    for i = 1:length(strs)
        str_out = strrep(str_out,strs{i},str);
    end
elseif ischar(strs)
    str_out = strrep(str_out,strs,str);
end
end