function rel_struct = replaceStringStruct(str1,str2,rel_struct)
rel_field = fieldnames(rel_struct);
for i = 1:length(rel_field)
    if ischar(rel_struct.(rel_field{i}))
        rel_struct.(rel_field{i}) = strrep(rel_struct.(rel_field{i}),str1,str2);
    elseif istable(rel_struct.(rel_field{i}))
        tab = rel_struct.(rel_field{i});
        for j = 1:size(tab,2)
            if ischar(tab{1,j})||(iscell(tab{1,j})&&ischar(tab{1,j}{:}))
                tab{:,j} = strrep(tab{:,j},str1,str2);
            end
        end
        rel_struct.(rel_field{i}) = tab;
    elseif isstruct(rel_struct.(rel_field{i}))
        rel_struct.(rel_field{i}) = replaceStringStruct(str1,str2,rel_struct.(rel_field{i}));
    end
end
end