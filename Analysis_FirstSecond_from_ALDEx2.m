 function [effect] = Analysis_FirstSecond_from_ALDEx2(folder_path,A,effect_threshold)

dirData = dir([folder_path '/ALDEx2_res' '/*_res.txt']);
All_file_name = {dirData.name};
for i = 1 : length(All_file_name)
    temp = regexp(All_file_name{i}, '_', 'split');
    if strcmp(temp{2},'Firstpart')
        
        [Diff_taxa,species_effect] = Read_ALDEx2_results([folder_path '/ALDEx2_res/' All_file_name{i}],effect_threshold);
        effect = zeros(1,size(A,1));
        [a,b,c] = intersect(1:size(A,1),species_effect(:,1),'stable');
        effect(b) = species_effect(:,2);
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