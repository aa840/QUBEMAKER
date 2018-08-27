function  script_getbondedparams( folder, N )

%BONDS
%Give the bonds within the molecule - changes should be made to all values

%Input
inputfolder  = horzcat('../Input_File', folder, '/');

%Output
fid = fopen(horzcat('../Output_File',folder,'bonds_ddec'), 'w');

%Get values necessary

charm_bonds_values = importdata(horzcat('../OPLS_Files/','charm_input_bonds'));
original_psf = importdata(horzcat(inputfolder,'original_AA_psf'));
new_psf = importdata(horzcat(inputfolder,'new_AA_psf'));
bonds = importdata(horzcat(inputfolder, 'bonds' )); 

%Get  OPLS number to name list and format correctly 
OPLS_number_to_name = importdata('../OPLS_Files/Number_to_Atom_type');
OPLS_number_to_name = strtrim(OPLS_number_to_name);

for i = 1:size(OPLS_number_to_name,1)
    OPLS_number_to_name{i} = strsplit(OPLS_number_to_name{i});
end

%Get Seminario values
seminario_bonds = importdata(horzcat('../Modified_Seminario_Files/','Average_Modified_Seminario_Bonds_All_AA'));
seminario_bonds_names = seminario_bonds.textdata;
seminario_bonds_values = seminario_bonds.data;
for i = 1:size(seminario_bonds_names,1)
    seminario_bonds_names{i} = strsplit(seminario_bonds_names{i}, '-');
    seminario_bonds_names{i} = strtrim(seminario_bonds_names{i});
end


bonds_together =  zeros(2,1);

%changes format
k = 1;
for i = 1:size(bonds,1)
    for j = 1:4
        if bonds(i, 2 *j ) > N || isnan(bonds(i, 2 *j ))
            break 
        end
        bonds_together(1,k) = bonds(i, 2 *j - 1);
        bonds_together(2,k) = bonds(i, 2 *j  );
        k = k + 1;
    end
end

number_bonds = size(bonds_together,2); 

%Get ids for the bonds listed for orig names
bond_id_orig = cell(number_bonds,2);
for i = 1:number_bonds
    bond_id_orig{i, 1} = original_psf.textdata{bonds_together(1,i),6};
    bond_id_orig{i, 2} = original_psf.textdata{bonds_together(2,i),6};
    
    %Make alphabetical 
    bond_id_orig(i, :) =  sort(bond_id_orig(i,:));
end

%Get ids for the bonds listed for new names
bond_id_new = cell(number_bonds,2);
for i = 1:number_bonds
    bond_id_new{i, 1} = new_psf.textdata{bonds_together(1,i),6};
    bond_id_new{i, 2} = new_psf.textdata{bonds_together(2,i),6};
    
    %Make alphabetical
    bond_id_orig(i, :) =  sort(bond_id_orig(i,:));
end

%Get OPLS atom types from CHARMM number format
OPLS_atom_names = cell(number_bonds,2);
for i = 1:number_bonds
    %Format of number and no letter
    OPLS_atom_names{i, 1} = bond_id_orig{i, 1}(2:4);
    OPLS_atom_names{i, 2} = bond_id_orig{i, 2}(2:4);
    
    j = 1;
    %Corresponding OPLS Names
    for j = 1:size(OPLS_number_to_name,1)
        if strcmp(char(OPLS_number_to_name{j}(1)), OPLS_atom_names{i, 1})
            OPLS_atom_names{i, 1}  = char(OPLS_number_to_name{j}(2));
        end
        if strcmp(char(OPLS_number_to_name{j}(1)), OPLS_atom_names{i, 2})
            OPLS_atom_names{i, 2}  = char(OPLS_number_to_name{j}(2));
        end
    end
    
end
%}

%Get bond parameters from this list - Original OPLS Values
bond_params= zeros(number_bonds,2);
n = 1; 
for i = 1:number_bonds
    for j = 1: size(charm_bonds_values.textdata,1)
        if strcmp(charm_bonds_values.textdata{j, 1},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_bonds_values.textdata{j, 2} , char( bond_id_orig(i, 2) ) )
            bond_params(n, : )  = charm_bonds_values.data(j, :);
            n = n + 1; 
        end
    end
end

number = []; 
%Add new parameters from Modified Seminario method
for i = 1:number_bonds
    for j = 1:size(seminario_bonds_names,1)
        if ( strcmp(OPLS_atom_names{i, 1}, seminario_bonds_names{j}(1)) && strcmp(OPLS_atom_names{i, 2}, seminario_bonds_names{j}(2)) ) || ( strcmp(OPLS_atom_names{i, 1}, seminario_bonds_names{j}(2)) && strcmp(OPLS_atom_names{i, 2}, seminario_bonds_names{j}(1)) )
             bond_params(i, :) = seminario_bonds_values(j,:);    
             number(i) = 1;
        end
    end
end

%Print to File
for i = 1:number_bonds
      fprintf(fid, '%s %s   %5.2f      %4.3f\n',  bond_id_new{i,1}, bond_id_new{i,2},  bond_params(i, 1 ) ,  bond_params(i, 2 ) );  
end

fclose(fid);

end
