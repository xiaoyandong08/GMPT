function Write_file_for_LEFSe(filename,OTU_table,group_name,subject_ID,taxa_name)


fib = fopen(filename,'wt');
fprintf(fib,'%s','OTU_ID');
for i = 1 : length(group_name)
    fprintf(fib,'%s\t','');
    fprintf(fib,'%s',group_name{i});
end
fprintf(fib,'%s\n','');
fprintf(fib,'%s','subject_ID');
for i = 1 : length(subject_ID)
    fprintf(fib,'%s\t','');
    fprintf(fib,'%s',subject_ID{i});
end
fprintf(fib,'%s\n','');

for i=1:size(OTU_table,1)
    fprintf(fib,'%s',num2str(taxa_name{i}));
    for j=1:size(OTU_table,2)
        fprintf(fib,'%s\t','');
        fprintf(fib,'%s',num2str(OTU_table(i,j)));
    end
    fprintf(fib,'%s\n','');

end
fclose(fib);

end