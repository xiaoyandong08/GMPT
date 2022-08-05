function [Diff_taxa,ALL_present,Fraction_present,effect] = Analysis_results_from_ANCOM(folder_path,N,part_pheno_group_name,group_name,wscore_threshold)

dirData = dir([folder_path '/ANCOM_res' '/*_res.txt']);

All_file_name = {dirData.name};
pair_name = nchoosek(part_pheno_group_name,2);
all_pair_name = nchoosek(group_name,2);
pair_index = [];
for i = 1 : size(pair_name,1)
    for j = 1 : size(all_pair_name,1)
        if (strcmp(pair_name{i,1},all_pair_name{j,1})&&strcmp(pair_name{i,2},all_pair_name{j,2}))
            pair_index = [pair_index j];
            break;
        end
    end
end

if length(pair_index) - size(pair_name,1)~=0
    disp('Wrong peojection of pair name')
end

count = 1;
for i = 1 : length(All_file_name)
    temp = regexp(All_file_name{i}, '_', 'split');
    index = str2num(temp{2});
    if strcmp(temp{2},'Firstpart')
        aa = 1;
    elseif find(pair_index==index)
        if index == 5
            aa = 1;
        end
        [Diff_taxa{index},effect{index}] = Read_ANCOM_results(N,[folder_path '/ANCOM_res/' All_file_name{i}],wscore_threshold);

    end
end


ALL_present = [];
Fraction_present = zeros(N,1);
for i = 1 : N
    count = 0;
    for j = 1 :  length(Diff_taxa)
        if length(intersect(i,Diff_taxa{1,j}(:)))>0
            count = count + 1;
        end
    end
    Fraction_present(i) = count;
    if count == length(Diff_taxa)
        ALL_present = [ALL_present i];
    end
end
end


