function [sum_charge_before, sum_charge_after] = process_ionised_psf( inputfolder,N )
% The psf file contains the bonds, angles, dihedrals, impropers and the old charges that need to be changed to ddec values.  
% This function creates a .psf file with new atom types and prints out the
% bonds, angles, impropers and dihedrals to seperate files

%Input
fileID = fopen(horzcat(inputfolder, 'ionized.psf'));

lj = importdata(horzcat(inputfolder,'lj_updated_eps80.dat'));

%Gets  bonds, angles, dihedrals
file_bonds = fopen(horzcat(inputfolder, 'bonds'), 'wt');
file_angles = fopen(horzcat(inputfolder, 'angles'), 'wt');
file_dihedrals = fopen(horzcat(inputfolder, 'dihedrals'), 'wt');
file_impropers = fopen(horzcat(inputfolder, 'impropers'), 'wt');

%The new ionized psf - changed charge and atom type 
file_new_ionized = fopen(horzcat(inputfolder,'new_ionized.psf'), 'wt');

%Just the amino acids part of the psf file
file_old_AA_psf = fopen(horzcat(inputfolder,'original_AA_psf'), 'wt');
file_new_AA_psf = fopen(horzcat(inputfolder,'new_AA_psf'), 'wt');

tline = fgets(fileID);

sum_charge_before = 0; 
sum_charge_after = 0; 

while ischar(tline)
    
    fprintf(file_new_ionized,'%s',tline);
    split_line = strsplit(tline); 

    if size(split_line,2) > 2 && strcmp(split_line{3}, '!NATOM')
        for i = 1:N
                %Number in atom type
                atomnumber = i - 1; 
                tline = fgets(fileID);
                split_line = strsplit(tline); 
                
                % Lengths fixed 
                if size(split_line{2},2) == 1
                    split_line{2} = horzcat(split_line{2}, '   ');
                end
                
                if size(split_line{2},2) == 2
                    split_line{2} = horzcat(split_line{2}, '  ');
                end
                
                if size(split_line{2},2) == 3
                    split_line{2} = horzcat(split_line{2}, ' ');
                end
                
                % Lengths fixed 
                if size(split_line{4},2) == 1
                    split_line{4} = horzcat(split_line{4}, ' ');
                end
                
                %Finds charge sum for OPLS
                sum_charge_before = sum_charge_before + str2num(split_line{8});
                %Replace charge
                split_line{8} = num2str(lj(i,1));
                
                %Shorten name
                tmp =  split_line{6}; 
                if  size(tmp,1) > 3
                    split_line{6} = tmp(1:3); 
                end
                if  size(tmp,1) == 2
                    split_line{6} = horzcat(tmp(1:2), ' ');
                end
                if  size(tmp,1) == 1
                    split_line{6} = horzcat(tmp(1), '  ');
                end
                
                %Finds charge sum for DDEC
                sum_charge_after = sum_charge_after + str2num(split_line{8});
                
                atomname = horzcat(split_line{7}(1), num2str(600 + atomnumber));
                
                if size(split_line{9}, 2) == 7
                    fprintf(file_new_ionized,'      %2s %s   %s    %s  %3s  %s  % 6.4f       % 6.4f           %s\n',split_line{2},split_line{3},split_line{4},split_line{5},split_line{6},atomname,str2num(split_line{8}), str2num(split_line{9}),split_line{10});
                else
                    fprintf(file_new_ionized,'      %2s %s   %s    %s  %3s  %s  % 6.4f        % 6.4f           %s\n',split_line{2},split_line{3},split_line{4},split_line{5},split_line{6},atomname,str2num(split_line{8}), str2num(split_line{9}),split_line{10});      
                end
                
                fprintf(file_new_AA_psf,'      %2s %s   %s    %s  %3s  %s  % f       %7s           %s\n',split_line{2},split_line{3},split_line{4},split_line{5},split_line{6},atomname,str2num(split_line{8}),split_line{9},split_line{10});
                fprintf(file_old_AA_psf,'%s',tline);
        end        
    end
    
    %Print Bonds to File
    if size(split_line,2) > 3 && strcmp(split_line{4}, 'bonds')
        tline = fgets(fileID);
        split_line = strsplit(tline); 
        while (str2num(split_line{2})  <= N )
              fprintf(file_new_ionized,'%s', tline);
              fprintf(file_bonds,'%s', tline);
              tline = fgets(fileID);
              split_line = strsplit(tline); 
        end
        fprintf(file_new_ionized,'%s', tline);
    end
    
    %Print Angles to File
    if size(split_line,2) > 3 && strcmp(split_line{4}, 'angles')
        tline = fgets(fileID);
        split_line = strsplit(tline); 
        while (str2num(split_line{2})  <= N )
              fprintf(file_new_ionized,'%s', tline);
              fprintf(file_angles,'%s', tline);
              tline = fgets(fileID);
              split_line = strsplit(tline); 
        end
        fprintf(file_new_ionized,'%s', tline);
    end
    
    %Print Dihedrals to File
    if size(split_line,2) > 3 && strcmp(split_line{4}, 'dihedrals')
        tline = fgets(fileID);
        split_line = strsplit(tline); 
        while (str2num(split_line{2})  <= N )
              fprintf(file_new_ionized,'%s', tline);
              fprintf(file_dihedrals,'%s', tline);
              tline = fgets(fileID);
              split_line = strsplit(tline); 
        end
        fprintf(file_new_ionized,'%s', tline);
    end
        
    %Print Impropers to File
    if size(split_line,2) > 3 && strcmp(split_line{4}, 'impropers')
        tline = fgets(fileID);
        split_line = strsplit(tline); 
        while (str2num(split_line{2})  <= N )
              fprintf(file_new_ionized,'%s', tline);
              fprintf(file_impropers,'%s', tline);
              tline = fgets(fileID);
              split_line = strsplit(tline); 
        end
        fprintf(file_new_ionized,'%s', tline);
    end
    
    tline = fgets(fileID);
    
end

end

