function [Diff_taxa,ALL_present,Fraction_present,effect] = Analysis_results_from_ALDEx2(folder_path,N,part_pheno_group_name,group_name,effect_threshold)

dirData = dir([folder_path '/ALDEx2_res' '/*_res.txt']);
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
%         [Diff_taxa{count},species_effect] = Read_ALDEx2_results([folder_path '/ALDEx2_res/' All_file_name{i}],effect_threshold);
%         effect{count} = zeros(1,N);
%         [a,b,c] = intersect(1:N,species_effect(:,1),'stable');
%         effect{count}(b) = species_effect(:,2);
%         count = count + 1;

        [Diff_taxa{index},species_effect] = Read_ALDEx2_results([folder_path '/ALDEx2_res/' All_file_name{i}],effect_threshold);
        effect{index} = zeros(1,N);
        [a,b,c] = intersect(1:N,species_effect(:,1),'stable');
        effect{index}(b) = species_effect(:,2);

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

function [Diff_taxa,effect] = Read_ALDEx2_results(filename,effect_threshold)
Diff_taxa = zeros(0,1);
fid = fopen(filename);
nextline = fgetl(fid);
nextline = fgetl(fid);

effect = [];
count = 1;
while isstr(nextline) & length(nextline)>0
    temp = regexp(nextline, '\t', 'split');
    temp2 = regexp(temp{1}, '"', 'split');
    index = str2num(temp2{2});
    effect = [effect; [index str2num(temp{15})]];
    if str2num(temp{15})>effect_threshold||str2num(temp{15})<-effect_threshold
        temp1 = regexp(temp{1}, '"', 'split');
        Diff_taxa(count,1) = str2num(temp1{2});
        count = count + 1;
    end
    nextline = fgetl(fid);
end
fclose(fid);

end
