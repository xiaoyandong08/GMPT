function Write_file_for_f1_macro_f1_micro(filename,data,threshold)
[taxa_size,~] = size(data);

fib = fopen(filename,'wt');
fprintf(fib,'%s\t%s\t%s\n','threshold','f1_macro','f1_micro');

for i=1:taxa_size

    fprintf(fib,'%s',[num2str(threshold(i))]);
    for j=1:size(data,2)
        fprintf(fib,'%s\t','');
        fprintf(fib,'%.6f',data(i,j));
    end
    fprintf(fib,'%s\n','');

end
fclose(fib);

end