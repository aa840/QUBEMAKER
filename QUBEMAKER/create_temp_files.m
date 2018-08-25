function create_temp_files(folder, N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Original atom types and bonds from PSF file

fid = fopen(horzcat(folder, 'ionized.psf'));
fidout = fopen(horzcat(folder, 'First_part_psf'), 'wt');
fidbond = fopen(horzcat(folder, 'Bonds_psf'), 'wt');
tline = fgets(fid);

while ischar(tline)
    if(strcmp(strtrim(tline(10:end)), strtrim('!NATOM')))
        tline = fgets(fid);
        for i = 1:N
            fprintf(fidout, '%s', tline);
            tline = fgets(fid);
        end
    end
    
    if(strcmp(strtrim(tline(10:end)), strtrim('!NBOND: bonds')))
        tline = fgets(fid);
        split = strsplit(tline);
        while size(tline,2) > 2 && (str2num(split{2}) <= N )
            fprintf(fidbond, '%s', tline);
            tline = fgets(fid);
            split = strsplit(tline);
        end
        
    end
    
    tline = fgets(fid);
end

fclose(fid);
fclose(fidout);

%Charges from .onetep file 

fid = fopen(horzcat(folder, 'ddec.onetep'));
fidout = fopen(horzcat(folder, 'charges.dat'), 'wt');
tline = fgets(fid);
while ischar(tline)
    if(strcmp(strtrim(tline), strtrim('DDEC Charges (X=0.02)')))
        tline = fgets(fid);
        tline = fgets(fid);
        tline = fgets(fid);
        tline = fgets(fid);
        for i = 1:N
            split_charges = strsplit(tline);
            fprintf(fidout, '%s   %s\n', split_charges{2}, split_charges{5});
            tline = fgets(fid);
        end
    end
   tline = fgets(fid);
end

fclose(fid); 
fclose(fidout); 

%AIM Volumes from .onetep file

fid = fopen(horzcat(folder, 'ddec.onetep'));
fidout = fopen(horzcat(folder, 'temp'), 'wt');
tline = fgets(fid);
while ischar(tline)
    if(strcmp(strtrim(tline), strtrim('DDEC Radial Moment Order  3 (X=0.02)')))
        tline = fgets(fid);
        tline = fgets(fid);
        tline = fgets(fid);
        tline = fgets(fid);
        for i = 1:N
            fprintf(fidout, '%s', tline);
            tline = fgets(fid);
        end
    end
   tline = fgets(fid);
end

end

