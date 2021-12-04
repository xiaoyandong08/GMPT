function Write_file_for_Spearman_correlation(filename,data,taxa_name)
[taxa_size,~] = size(data);

fib = fopen(filename,'wt');
fprintf(fib,'%s\t%s\t%s\n','OTU_ID','Spearman','Interaction');

for i=1:taxa_size
    if taxa_name(i) == 1
        fprintf(fib,'%s',['taxa' 'X']);
    else
        fprintf(fib,'%s',['taxa' num2str(taxa_name(i))]);
    end
    for j=2:size(data,2)
        fprintf(fib,'%s\t','');
        fprintf(fib,'%.6f',data(i,j));
    end
    fprintf(fib,'%s\n','');

end
fclose(fib);

end