function Write_file_for_OTU_table(filename,OTU_table,group_name)
[taxa_size,sample_size] = size(OTU_table{1});
OTU_table = cell2mat(OTU_table);

fib = fopen(filename,'wt');
fprintf(fib,'%s','OTU_ID');
for i = 1 : length(group_name)
    fprintf(fib,'%s\t','');
    for j = 1 : sample_size-1
        fprintf(fib,'%s\t',[group_name{i} '_sample' num2str(j)]);
    end
    fprintf(fib,'%s',[group_name{i} '_sample' num2str(j+1)]);
end
fprintf(fib,'%s\n','');

for i=1:taxa_size
    if i == 1
        fprintf(fib,'%s',['taxa' 'X']);
    else
        fprintf(fib,'%s',['taxa' num2str(i)]);
    end
    for j=1:size(OTU_table,2)
        fprintf(fib,'%s\t','');
        fprintf(fib,'%.6f',OTU_table(i,j));
    end
    fprintf(fib,'%s\n','');

end
fclose(fib);

end