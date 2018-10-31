function get_mass( folder, N )
%Extracts the masses and prints to file

%Input
inputfolder  = horzcat('../Input_File', folder ,  '/' );

%Output
fid = fopen(horzcat('../Output_File',folder,'masses'), 'w');

%Get values necessary
new_psf = importdata(horzcat(inputfolder,'new_AA_psf'));

masses = cell(1,N);

fprintf(fid, '%s \n', 'ATOMS');

%Get ids for the bonds listed for new names
for i = 1:N %All angles
    names = new_psf.textdata{i, 6};
    switch names(1)
        case 'N'
            masses{i} = '14.007000';
        case 'S'
            masses{i} = '32.060000';
        case 'O'
            masses{i} = '15.999400';
        case 'C'
            masses{i} = '12.011000';
        case 'H'
            masses{i} = '1.008000';
    end
    fprintf(fid, 'MASS %d %s  %s \n',  i, new_psf.textdata{i, 6}, masses{i});
end   
    
fprintf(fid, '%s %d %s \n', 'MASS', (N + 1) ,'HT   1.008000');
fprintf(fid, '%s %d %s \n', 'MASS', (N + 2) ,'OT   15.999400');
fprintf(fid, '%s %d %s \n', 'MASS', (N + 3) ,'SOD   22.989800');
fprintf(fid, '%s %d %s \n', 'MASS', (N + 4) ,'CLA  35.450000');

fclose(fid);   

end

