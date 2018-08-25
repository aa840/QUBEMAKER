function  script_angleparams( folder,  N )

%ANGLES
%Give the angles within the molecule - changes should be made to all values

%Input
inputfolder  = horzcat('../Input_File', folder ,  '/' );

%Output
fid = fopen(horzcat('../Output_File',folder,'angles_ddec'), 'w');

%Get values necessary
charm_angles_values = importdata(horzcat('../OPLS_Files/','charm_input_angles'));
original_psf = importdata(horzcat(inputfolder,'original_AA_psf'));
new_psf = importdata(horzcat(inputfolder,'new_AA_psf'));
angle = importdata(horzcat(inputfolder, 'angles' ));


%Get  OPLS number to name list and format correctly
OPLS_number_to_name = importdata('../OPLS_Files/Number_to_Atom_type');
OPLS_number_to_name = strtrim(OPLS_number_to_name);

for i = 1:size(OPLS_number_to_name,1)
    OPLS_number_to_name{i} = strsplit(OPLS_number_to_name{i});
end

%Get Seminario values
seminario_angles = importdata(horzcat('../Modified_Seminario_Files/','Average_Modified_Seminario_Angles_All_AA'));
seminario_angles_names = seminario_angles.textdata;
seminario_angles_value = seminario_angles.data;

for i = 1:size(seminario_angles_names,1)
    seminario_angles_names{i} = strsplit(seminario_angles_names{i}, '-');
    seminario_angles_names{i} = strtrim(seminario_angles_names{i});
end


angles_together =  zeros(3,1);

%Changes format
k = 1;
for i = 1:size(angle,1)
    for j = 1:3
        if angle(i, 3 *j ) > N || isnan(angle(i, 3 *j ))
            break
        end
        angles_together(1,k) = angle(i, 3 *j - 2);
        angles_together(2,k) = angle(i, 3 *j  - 1);
        angles_together(3,k) = angle(i, 3 *j );
        k = k + 1;
    end
end

number_angle = size(angles_together,2);

%Get ids for the bonds listed for orig names
angle_id_orig = cell(number_angle,3);
for i = 1:number_angle %All angles
    angle_id_orig{i, 1} = original_psf.textdata{angles_together(1,i),6};
    angle_id_orig{i, 2} = original_psf.textdata{angles_together(2,i),6};
    angle_id_orig{i, 3} = original_psf.textdata{angles_together(3,i),6};
end

%Get ids for the bonds listed for new names
angle_id_new = cell(number_angle,3);
for i = 1:number_angle %All angles
    angle_id_new{i, 1} = new_psf.textdata{angles_together(1,i),6};
    angle_id_new{i, 2} = new_psf.textdata{angles_together(2,i),6};
    angle_id_new{i, 3} = new_psf.textdata{angles_together(3,i),6};
end

%Get OPLS atom types from CHARMM number format
OPLS_atom_names = cell(number_angle,3);
for i = 1:number_angle
    %Format of number and no letter
    OPLS_atom_names{i, 1} = angle_id_orig{i, 1}(2:4);
    OPLS_atom_names{i, 2} = angle_id_orig{i, 2}(2:4);
    OPLS_atom_names{i, 3} = angle_id_orig{i, 3}(2:4);
    
    %Corresponding OPLS Names
    j = 1; 
    changed = 'n';
    while changed == 'n'
        if strcmp(char(OPLS_number_to_name{j}(1)), OPLS_atom_names{i, 1})
            OPLS_atom_names{i, 1}  = char(OPLS_number_to_name{j}(2));
            changed = 'y';
        end
        if strcmp(char(OPLS_number_to_name{j}(1)), OPLS_atom_names{i, 2})
            OPLS_atom_names{i, 2}  = char(OPLS_number_to_name{j}(2));
            changed = 'y';
        end
        if strcmp(char(OPLS_number_to_name{j}(1)), OPLS_atom_names{i, 3})
            OPLS_atom_names{i, 3}  = char(OPLS_number_to_name{j}(2));
            changed = 'y';
        end
        j = j + 1; 
    end
end

%Get angle parameters from this list - Original OPLS Values
angle_params= zeros(number_angle,2);
n = 1;
for i = 1:size(angles_together, 2) % all angles
    for j = 1:size(charm_angles_values.textdata,1) %all charm entries
        if ( strcmp(charm_angles_values.textdata{j, 1},  char(angle_id_orig(i, 1) ) ) && strcmp( charm_angles_values.textdata{j, 2} , char( angle_id_orig(i, 2) ) ) && strcmp( charm_angles_values.textdata{j, 3} , char( angle_id_orig(i, 3) ) ) )  || ...
                ( strcmp(charm_angles_values.textdata{j, 3},  char(angle_id_orig(i, 1) ) ) && strcmp( charm_angles_values.textdata{j, 2} , char( angle_id_orig(i, 2) ) ) && strcmp( charm_angles_values.textdata{j, 1} , char( angle_id_orig(i, 3) ) ) )
            angle_params(n, : )  = charm_angles_values.data(j, :);
            n = n + 1;
        end
    end
end

number = [];
%Add new parameters from Modified Seminario 
for i = 1:number_angle
    for j = 1:size(seminario_angles_value,1)
        if ( strcmp(OPLS_atom_names{i, 1}, seminario_angles_names{j}(1)) && strcmp(OPLS_atom_names{i, 2}, seminario_angles_names{j}(2)) && strcmp(OPLS_atom_names{i, 3}, seminario_angles_names{j}(3)) ) || ( strcmp(OPLS_atom_names{i, 1}, seminario_angles_names{j}(3)) && strcmp(OPLS_atom_names{i, 2}, seminario_angles_names{j}(2)) && strcmp(OPLS_atom_names{i, 3}, seminario_angles_names{j}(1)) )
            angle_params(i, :) = seminario_angles_value(j,:);
            number(i) = 1;
        end
    end
end

%Print to File
for i = 1:number_angle
    if angle_params(i,1) > 100
        fprintf(fid, '%s %s %s  % 4.2f   %5.2f \n',  angle_id_new{i,1}, angle_id_new{i,2}, angle_id_new{i,3}, angle_params(i,1) , angle_params(i,2)  );
    else     
        fprintf(fid, '%s %s %s   % 4.2f   %5.2f \n',  angle_id_new{i,1}, angle_id_new{i,2}, angle_id_new{i,3}, angle_params(i,1) , angle_params(i,2)  );
    end
end

fclose(fid);

end

