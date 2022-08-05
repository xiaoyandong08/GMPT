 function [effect] = Analysis_FirstSecond_from_ANCOM(folder_path,A,wscore_threshold)
N = size(A,1);
dirData = dir([folder_path '/ANCOM_res' '/*_res.txt']);
All_file_name = {dirData.name};
for i = 1 : length(All_file_name)
    temp = regexp(All_file_name{i}, '_', 'split');
    if strcmp(temp{2},'Firstpart')
        
        [Diff_taxa,effect] = Read_ANCOM_results(N,[folder_path '/ANCOM_res/' All_file_name{i}],wscore_threshold);

    end

    
end
 end
