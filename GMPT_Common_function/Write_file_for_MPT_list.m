function Write_file_for_MPT_list(filename,plot_matrix,phenotype_name,taxa_name)
[taxa_size,phenotype_size] = size(plot_matrix);

fib = fopen(filename,'wt');
fprintf(fib,'%s','OTU_ID');
for i = 1 : length(phenotype_name)
    fprintf(fib,'%s\t','');
    fprintf(fib,'%s',phenotype_name{i});
end
fprintf(fib,'%s\n','');

for i=1:taxa_size
    if taxa_name(i) == 1
        fprintf(fib,'%s',['taxa' 'X']);
    else
        fprintf(fib,'%s',['taxa' num2str(taxa_name(i))]);
    end
    for j=1:size(plot_matrix,2)
        fprintf(fib,'%s\t','');
        fprintf(fib,'%.6f',plot_matrix(i,j));
    end
    fprintf(fib,'%s\n','');

end
fclose(fib);

end