%% Convert cell to data structure


function data_struct = cell_to_struct(data,h,temp_fs,pre_post)

    data_struct = struct();
    data_struct.fsample = temp_fs;
    data_struct.time = repmat({linspace(pre_post(1),pre_post(end),size(data{1},2))},1,length(data));
    data_struct.label = arrayfun(@(x) ['chan',num2str(x+h-1)],1:size(data{1},1),'UniformOutput',false);

    data_struct.trial = data;
end