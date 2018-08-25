fid = fopen('Ordered', 'wt')

i = 1; 
j = 1; 
while j <= size(array,2)
    j
    i
    if  array(j) ~= 0
        fprintf(fid, '%s %s \n', num2str(j), number_atom_type{array(j), 2});
    else
        fprintf(fid, '%s %s \n', num2str(j), 'ER');        
    end
    
    j =j + 1; 
end

fclose(fid)