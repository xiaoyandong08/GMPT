function [Diff_taxa,effect] = Read_ANCOM_results(N,filename,wscore_threshold)
Diff_taxa = zeros(0,1);
fid = fopen(filename);
nextline = fgetl(fid);

temp = regexp(nextline, '\t', 'split');
for i = 1 : length(temp)
    temp1 = regexp(temp{i}, '"', 'split');
    temp1(cellfun(@isempty,temp1)) = [];
    if ~isempty(temp1)
        temp1 = regexp(temp1{1}, 'detected_', 'split');
        temp1(cellfun(@isempty,temp1)) = [];
        if str2num(temp1{1})==wscore_threshold
            threshold_index = i;
            break;
        end
    end
end


nextline = fgetl(fid);

effect = [];
count = 1;
while isstr(nextline) & length(nextline)>0
    temp = regexp(nextline, '\t', 'split');
    
    temp2 = regexp(temp{2}, '"', 'split');
    temp3 = regexp(temp2{2}, 'OTU', 'split');
    index = str2num(temp3{2});
    
    effect = [effect; [index str2num(temp{3}) 0]];
    if strcmp(temp{threshold_index},'TRUE')
        effect(end,3) = 1;
        Diff_taxa(count,1) = index;
        count = count + 1;
    end
    nextline = fgetl(fid);
end
fclose(fid);

Diff_taxa = sort(Diff_taxa,'ascend');
[~,index] = sort(effect(:,1),'ascend');
effect = effect(index,:);
if size(effect,1)~=N
    index = setdiff([1:N],effect(:,1));
    effect = [effect;[index' zeros(length(index),2)]];
    [~,index] = sort(effect(:,1),'ascend');
    effect = effect(index,:);
end
if sum(find(effect(:,3))-Diff_taxa)~=0
    error('Wrong!')
end
end