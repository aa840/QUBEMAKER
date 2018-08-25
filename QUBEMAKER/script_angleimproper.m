function  script_angleimproper( folder, N )

%Impropers
%Give the impropers within the molecule - OPLS values

%Input
inputfolder  = horzcat('../Input_File', folder , '/');

%Output
fid_out = fopen(horzcat('../Output_File',folder,'improper_ddec'), 'w');

%Log File
fid_log = fopen(horzcat('../Output_File/log_dihedrals'), 'wt');

%Get values necessary
charm_improper_values = importdata(horzcat('../OPLS_Files/','charm_input_improper'));
 
original_psf = importdata(horzcat(inputfolder,'original_AA_psf'));
new_psf = importdata(horzcat(inputfolder,'new_AA_psf'));
improper = importdata(horzcat(inputfolder, 'impropers' )); 
improper_together =  zeros(4,1);

%changes format to one improper per row
k = 1;
for i = 1:size(improper,1)
    for j = 1:2
        %Only want impropers in AA 
        if improper(i, 4 *j ) > N || isnan(improper(i, 4 *j ))
            break
        end
        improper_together(1,k) = improper(i, 4 *j - 3);
        improper_together(2,k) = improper(i, 4 *j  - 2);
        improper_together(3,k) = improper(i, 4 *j  - 1);
        improper_together(4,k) = improper(i, 4 *j  );
        k = k + 1;
    end
end

number_improper = size(improper_together,2); 

%get ids for the bonds listed for orig names
bond_id_orig = cell(number_improper,4);
for i = 1:number_improper %all improper
    bond_id_orig{i, 1} = original_psf.textdata{improper_together(1,i),6};
    bond_id_orig{i, 2} = original_psf.textdata{improper_together(2,i),6};
    bond_id_orig{i, 3} = original_psf.textdata{improper_together(3,i),6};
    bond_id_orig{i, 4} = original_psf.textdata{improper_together(4,i),6};
end

%get ids for the bonds listed for new names
bond_id_new = cell(number_improper,4);
for i = 1:number_improper %all impropers
    bond_id_new{i, 1} = new_psf.textdata{improper_together(1,i),6};
    bond_id_new{i, 2} = new_psf.textdata{improper_together(2,i),6};
    bond_id_new{i, 3} = new_psf.textdata{improper_together(3,i),6};
    bond_id_new{i, 4} = new_psf.textdata{improper_together(4,i),6};
    %make alphabetical
    %  bond_id_orig(i, :) =  sort(bond_id_orig(i,:));
end

%get dihedral parameters from this list
improper_params= zeros(number_improper,3);
n = 1;

array_i =  zeros(number_improper,1); %Number of terms in the improper

for i = 1:number_improper % all impropers
    no = 0;
    j = 1;
    n_before = n;
    while j <= size( charm_improper_values.textdata, 1) %all charm entries
        if ( strcmp(charm_improper_values.textdata{j, 1},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_improper_values.textdata{j, 2} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_improper_values.textdata{j, 3} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_improper_values.textdata{j, 4} , char( bond_id_orig(i, 4) ) ) || ...
                strcmp(charm_improper_values.textdata{j, 4},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_improper_values.textdata{j, 3} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_improper_values.textdata{j, 2} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_improper_values.textdata{j, 1} , char( bond_id_orig(i, 4) ) )) % or reverse order
            improper_params(n, : )  = charm_improper_values.data(j, :);
            n = n + 1;
            no = no + 1;
            j = j + 1;
            if ( strcmp(charm_improper_values.textdata{j, 1},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_improper_values.textdata{j, 2} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_improper_values.textdata{j, 3} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_improper_values.textdata{j, 4} , char( bond_id_orig(i, 4) ) ) || ...
                    strcmp(charm_improper_values.textdata{j, 4},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_improper_values.textdata{j, 3} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_improper_values.textdata{j, 2} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_improper_values.textdata{j, 1} , char( bond_id_orig(i, 4) ) )) % or reverse order
                improper_params(n, : )  = charm_improper_values.data(j , :);
                n = n + 1;
                no = no + 1;
                j = j + 1;
                if ( strcmp(charm_improper_values.textdata{j, 1},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_improper_values.textdata{j, 2} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_improper_values.textdata{j, 3} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_improper_values.textdata{j, 4} , char( bond_id_orig(i, 4) ) ) || ...
                        strcmp(charm_improper_values.textdata{j, 4},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_improper_values.textdata{j, 3} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_improper_values.textdata{j, 2} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_improper_values.textdata{j, 1} , char( bond_id_orig(i, 4) ) )) % or reverse order
                    improper_params(n, : )  = charm_improper_values.data(j, :);
                    n = n + 1;
                    no = no + 1;
                    j = j + 1;
                end
            end
        else
            no = 1;
            j = j + 1;
        end
    end
    
    %Below is case if improper not found
    if n_before == n
        improper_params(n, : ) = [0,0,0];
        n = n +  1;
        no = 1;
         fprintf(fid_log, 'PROBLEM with IMPROPER TORSIONALS  \n' );
    end

    array_i(i) = no;
end

m = 1;
for i = 1:size(array_i, 1)
    for j = 1:array_i(i) %for all improper lines
        if improper_params(m,1) <  10 
            fprintf(fid_out, '%s %s %s %s    % 6.4f   %1.0f  % 4.1f\n',  bond_id_new{i,1}, bond_id_new{i,2}, bond_id_new{i,3},  bond_id_new{i,4}, improper_params(m,1) , improper_params(m,2),  improper_params(m,3)  );
        else
            fprintf(fid_out, '%s %s %s %s   % 6.4f   %1.0f  % 4.1f\n',  bond_id_new{i,1}, bond_id_new{i,2}, bond_id_new{i,3},  bond_id_new{i,4}, improper_params(m,1) , improper_params(m,2),  improper_params(m,3)  );
        end
        m = m+1;
    end
end

fclose(fid_out);
fclose(fid_log);

end
