function script_angledihedral( folder, N, tp_name)

%DIHEDRAL

number_new_tp = 0; %Checking mechanism

%Input folders
inputfolder  = horzcat('./Input_File', folder , '/');

%Output printed
fid_out = fopen(horzcat('./Output_File',folder,'dihedral_ddec_', tp_name), 'wt');

%Log file
fid_log = fopen(horzcat('./Output_File/log_dihedrals'), 'wt');

%Get parameters necessary
charm_dihedral_values = importdata(horzcat('./OPLS_Files/','charm_input_dihedral'));
charm_dihedral_values_wildcards = importdata(horzcat('./OPLS_Files/','charm_input_dihedral_wildcards'));
copy_charm_dihedral_values = charm_dihedral_values;
original_psf = importdata(horzcat(inputfolder,'original_AA_psf'));
new_psf = importdata(horzcat(inputfolder,'new_AA_psf'));
dihedral = importdata(horzcat(inputfolder, 'dihedrals' ));

%Get backbone atom names
tmp = importdata('./New_torsional_parameters/name_torsion_params' );
name_backbone_tp = cell(size(tmp,1),4);

for i=1:size(tmp,1)
    name_backbone_tp(i,:)  = strsplit(strtrim(tmp{i}));
end

old_torsional_params = importdata('./OPLS_Files/OPLS_psi_psi_torsion_params' );
old_torsional_params  = old_torsional_params./2; %OPLS to CHARMM

%Backbone Terms
backbone_tp = importdata(horzcat('./New_torsional_parameters/', tp_name ));
backbone_tp = backbone_tp./2; %Change from OPLS to CHARMM

%Changes format to one dihedral per row
dihedral_together =  zeros(4,1);
k = 1;
for i = 1:size(dihedral,1)
    for j = 1:2

        if dihedral(i, 4 *j ) > N || isnan(dihedral(i, 4 *j ))
            break
        end
        
        dihedral_together(1,k) = dihedral(i, 4 *j - 3);
        dihedral_together(2,k) = dihedral(i, 4 *j  - 2);
        dihedral_together(3,k) = dihedral(i, 4 *j  - 1);
        dihedral_together(4,k) = dihedral(i, 4 *j  );
        k = k + 1;
    end
end

number_dihedrals = size(dihedral_together,2);
dihedrals_change = [];

numbers_charmm  = importdata('./OPLS_Files/Number_to_Atom_type');

number_atom_type = cell(size(numbers_charmm,1),2);

for i = 1:size(number_atom_type,1)
    number_atom_type(i, :) = strsplit(strtrim(numbers_charmm{i}));
end

%Get ids for the bonds listed for original names
bond_id_orig = cell(number_dihedrals,4);
for i = 1:number_dihedrals %all dihedral
    bond_id_orig{i, 1} = original_psf.textdata{dihedral_together(1,i),6};
    bond_id_orig{i, 2} = original_psf.textdata{dihedral_together(2,i),6};
    bond_id_orig{i, 3} = original_psf.textdata{dihedral_together(3,i),6};
    bond_id_orig{i, 4} = original_psf.textdata{dihedral_together(4,i),6};
end

%Get ids for the bonds listed for new names
bond_id_new = cell(number_dihedrals,4);
residue_id = cell(number_dihedrals,1);

%prev_residue = original_psf.textdata{dihedral_together(1,1),4};
residue_number_array = zeros(1,number_dihedrals);

%New residue names and residue types
for i = 1:number_dihedrals %all dihedral
    bond_id_new{i, 1} = new_psf.textdata{dihedral_together(1,i),6};
    bond_id_new{i, 2} = new_psf.textdata{dihedral_together(2,i),6};
    bond_id_new{i, 3} = new_psf.textdata{dihedral_together(3,i),6};
    bond_id_new{i, 4} = new_psf.textdata{dihedral_together(4,i),6};
    
    residue_id{i, 1} = original_psf.textdata{dihedral_together(1,i),4};
    residue_id{i, 2} = original_psf.textdata{dihedral_together(2,i),4};
    residue_id{i, 3} = original_psf.textdata{dihedral_together(3,i),4};
    residue_id{i, 4} = original_psf.textdata{dihedral_together(4,i),4};
    
    %Number of the residue 
    residue_number_array(i) = str2num(original_psf.textdata{dihedral_together(1,i),3});
end

%Total number of residues
residue_count = residue_number_array(end);


%Get dihedral parameters from this list
dihedral_params= zeros(number_dihedrals,3);
dihedral_params_index= zeros(number_dihedrals,1);

copy_dihedral_params= zeros(number_dihedrals,3);

fprintf(fid_log, '%s %s \n', 'Start log', folder);
fprintf(fid_log, '%s %s\n', 'Time is now: ', datestr(clock, 0));

n = 1;
array_no =  zeros(number_dihedrals,1); %The number of terms in the torsional parameter

c_alpha = cell(1,residue_count);
c_beta = cell(1,residue_count);

%Find old values and C_beta
for i = 1:number_dihedrals % all dihedral
    
    %Copied as changed when wildcards considered
    charm_dihedral_values = copy_charm_dihedral_values;
    
    no = 0;
    j = 1;
    
    n_before = n;
    
    changed = 'n';
    
    %Go through all charm entries to find
    while (changed == 'n') && j <= size(charm_dihedral_values.data, 1) && strcmp(changed, 'n')%all charm entries
        
        %Finds the  dihedral terms
        if ( strcmp(charm_dihedral_values.textdata{j, 1},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_dihedral_values.textdata{j, 2} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_dihedral_values.textdata{j, 3} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_dihedral_values.textdata{j, 4} , char( bond_id_orig(i, 4) ) ) || ...
                strcmp(charm_dihedral_values.textdata{j, 4},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_dihedral_values.textdata{j, 3} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_dihedral_values.textdata{j, 2} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_dihedral_values.textdata{j, 1} , char( bond_id_orig(i, 4) ) )) % or reverse order
            dihedral_params(n, : )  = charm_dihedral_values.data(j, :);
            dihedral_params_index(n) = j;
            copy_dihedral_params(n, : ) = charm_dihedral_values.data(j, :);
            
            changed = 'y';
            
            n = n + 1;
            no = no  + 1;
            
            j = j + 1;
            
            if ( strcmp(charm_dihedral_values.textdata{j, 1},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_dihedral_values.textdata{j, 2} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_dihedral_values.textdata{j, 3} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_dihedral_values.textdata{j, 4} , char( bond_id_orig(i, 4) ) ) || ...
                    strcmp(charm_dihedral_values.textdata{j, 4},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_dihedral_values.textdata{j, 3} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_dihedral_values.textdata{j, 2} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_dihedral_values.textdata{j, 1} , char( bond_id_orig(i, 4) ) )) % or reverse order
                dihedral_params(n, : )  = charm_dihedral_values.data(j , :);
                copy_dihedral_params(n, : ) = charm_dihedral_values.data(j, :);
                
                n = n + 1;
                no = no  + 1;
                j = j + 1;
                
                if (  strcmp(charm_dihedral_values.textdata{j, 1},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_dihedral_values.textdata{j, 2} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_dihedral_values.textdata{j, 3} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_dihedral_values.textdata{j, 4} , char( bond_id_orig(i, 4) ) ) || ...
                        strcmp(charm_dihedral_values.textdata{j, 4},  char(bond_id_orig(i, 1) ) ) && strcmp( charm_dihedral_values.textdata{j, 3} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_dihedral_values.textdata{j, 2} , char( bond_id_orig(i, 3) ) ) && strcmp( charm_dihedral_values.textdata{j, 1} , char( bond_id_orig(i, 4) ) )) % or reverse order
                    dihedral_params(n, : )  = charm_dihedral_values.data(j, :);
                    copy_dihedral_params(n, : ) = charm_dihedral_values.data(j, :);
                    
                    n = n + 1;
                    no = no  + 1;
                    j = j + 1;
                end
            end
        else
            j = j + 1;
        end
    end
    
    n = n - (no + 1);
    
    %Using the psi' torsional parameters find the c_beta and c_atom atoms 
    if ((n+3) <= size(dihedral_params,1)) && (dihedral_params((n + 1), 1) == old_torsional_params(4,1) ) && (dihedral_params((n + 2), 1) == old_torsional_params(4,2)  ) && (dihedral_params((n + 3), 1) == old_torsional_params(4,3) )
        tmp = bond_id_orig{i,1};

        if strcmp(tmp(1), 'C')
            c_beta{residue_number_array(i)} = {bond_id_new{i,1}};
            c_alpha{residue_number_array(i)} = {bond_id_new{i,2}};
        else
            c_beta{residue_number_array(i)} = {bond_id_new{i,4}};
            c_alpha{residue_number_array(i)} = {bond_id_new{i,3}};
        end
        
    end
    
    n = n + (no + 1);
    
    
    if n_before == n
        fprintf(fid_log, '%s %s\n', 'No term found for residue number:',  int2str(residue_number_array(i)) );

        for l = 1:4
            fprintf(fid_log, '%s \n', bond_id_orig{i, l});
        end
        
        j = 1;

        %Records if the wildcard is changed
        wc_changed = 'n';
        while j <= size(charm_dihedral_values_wildcards.data, 1) &&  strcmp(wc_changed, 'n')%all charm entries
            
            %Finds the  dihedral terms
            if ( strcmp( charm_dihedral_values_wildcards.textdata{j, 2} , char( bond_id_orig(i, 2) ) ) && strcmp( charm_dihedral_values_wildcards.textdata{j, 3} , char( bond_id_orig(i, 3) ) ) )
                dihedral_params(n, : )  = charm_dihedral_values_wildcards.data(j, :);
                copy_dihedral_params(n, : ) = charm_dihedral_values_wildcards.data(j, :);
                
                wc_changed = 'y';
                
                n = n + 1;
                no = no  + 1;
                
                j = j + 1;
                
                        fprintf(fid_log, '%s\n', 'Wildcards term found!');
            else
                j = j + 1;
            end
        end
        if strcmp(wc_changed,'n')
            fprintf(fid_log, '%s\n', ' NO WILDCARDS FOUND - PROBLEM!');
            dihedral_params(n, : ) = [0,0,0];
            n = n + 1;
            no = 1;
        end
        
    end
    
    array_no(i) = no;
    
end

n=1;

%Finding new values
new_dihedral_parameter = [];

replaced = 'no';

new_array_no = array_no;

%Replace terms

for i = 1:number_dihedrals % all dihedral
    
    %Number of terms in dihedral parameters
    no = array_no(i);
    
    %Change TP to new values
    [ OPLS_names ] = opls_atom_type_to_charmm( {bond_id_orig{i,:}}, number_atom_type );

    %End of protein treated like rest of terms
    for ii = 1:4
        if strcmp(OPLS_names{ii}, 'N3')
            OPLS_names{ii} = 'N';
        end
    end
     
    j = 1;

    %Take the middle residue value
    residue_name = residue_id{i,2};
    
    %Get residue torsional parameters for current residue
    %if Proline need to replace psi psi' too
    switch residue_name
        case {'SER','THR', 'PRO', 'ALA', 'GLY'}
            sidechains_values = importdata(horzcat('./New_torsional_parameters/sidechains_values/', residue_id{i,2} ));
            sidechains_values = sidechains_values./2; %Change from OPLS to Torsional
            
            tmp = importdata(horzcat( './New_torsional_parameters/sidechains_names/', residue_id{i,2} ));
            sidechains_replace = cell(size(tmp,1),4);

            for m=1:size(tmp,1)
                sidechains_replace(m,:)  = strsplit(strtrim(tmp{m}));
            end
            
            %All torsional parameters are included in the sidechain file
            %for proline, serine and threonine, alanine and glycine
            values_tp_replace = sidechains_values;
            names_tp_replace= sidechains_replace;
            
        otherwise
            sidechains_values = importdata(horzcat( './New_torsional_parameters/sidechains_values/', residue_id{i,2} ));
            sidechains_values = sidechains_values./2; %Change from OPLS to Torsional
            
            tmp = importdata(horzcat('./New_torsional_parameters/sidechains_names/', residue_id{i,2} ));
            sidechains_replace = cell(size(tmp,1),4);

            for m=1:size(tmp,1)
                sidechains_replace(m,:)  = strsplit(strtrim(tmp{m}));
            end
            
            values_tp_replace = [backbone_tp; sidechains_values];
            names_tp_replace= [ name_backbone_tp; sidechains_replace];
            
    end
    
    %Loop through all new torsional parameters and replace as necessary
    while j <= size(values_tp_replace,1)

        %If one of the values to replace then replaced
        if strcmp(OPLS_names{1}, names_tp_replace{j,1}) && strcmp(OPLS_names{2}, names_tp_replace{j,2}) && strcmp(OPLS_names{3}, names_tp_replace{j,3}) && strcmp(OPLS_names{4}, names_tp_replace{j,4})
            n = n - 1; %Start beginning dihedral parameters
            
            %Backbone test - Must have OPLS-AA/M's params as old
            if ( sum(strcmp( residue_id{i,1}, {'SER','THR', 'PRO', 'ALA', 'GLY'})) || (j > 4) ||( j <=4 && (dihedral_params((n + 1), 1) == old_torsional_params(j,1) ) && (dihedral_params((n + 2), 1) == old_torsional_params(j,2)  ) && (dihedral_params((n + 3), 1) == old_torsional_params(j,3) ))  )
                
                %Sidechain test - Must have C-beta term for one atom
                if ( j <=4 || (residue_number_array(i) == residue_count) || ( j > 4 && (strcmp({bond_id_new{i,1}}, c_beta{residue_number_array(i)}) || strcmp({bond_id_new{i,2}}, c_beta{residue_number_array(i)}) || strcmp({bond_id_new{i,3}}, c_beta{residue_number_array(i)})  || strcmp({bond_id_new{i,4}}, c_beta{residue_number_array(i)}) ) ))
                    m = size(new_dihedral_parameter,1);

                    for k = 1:3 %Always have 3 terms
                        new_dihedral_parameter(m+k,1) = values_tp_replace(j,k);
                        new_dihedral_parameter(m+k,2)  =  k;
                        
                        replaced = 'yes';
                        new_array_no(i) = 3;
                    end
                    
                    %Phase of cosine
                    new_dihedral_parameter(m + 1,3) = 0;
                    new_dihedral_parameter(m + 2,3) = 180;
                    new_dihedral_parameter(m + 3,3) = 0;
                    
                    number_new_tp = number_new_tp + 1;
                    
                    %If replaced then move on to next dihedral value
                    j = size(values_tp_replace,1) + 1;
                end
            end
            n = n +  1;
        end
        
        %Backwards directions
        if j < size(names_tp_replace,1) && strcmp(OPLS_names{1}, names_tp_replace{j,4}) && strcmp(OPLS_names{2}, names_tp_replace{j,3}) && strcmp(OPLS_names{3}, names_tp_replace{j,2}) && strcmp(OPLS_names{4}, names_tp_replace{j,1})
            
            n= n - 1; %Start beginning dihedral parameters
            
            %Tests if proline, serine, threonine, a sidechain, or a backbone
            if ( sum(strcmp( residue_id{i,1}, {'SER','THR', 'PRO', 'ALA', 'GLY'})) || (j > 4) ||( j <= 4 && (dihedral_params((n + 1), 1) == old_torsional_params(j,1) ) && (dihedral_params((n + 2), 1) == old_torsional_params(j,2)  ) && (dihedral_params((n + 3), 1) == old_torsional_params(j,3) ))  )
                
                %Sidechain test - Must have C-beta term in unless at end as
                %not identified with current method
                if ( j <=4 || (residue_number_array(i) == residue_count) || ( j > 4 && (strcmp({bond_id_new{i,1}}, c_beta{residue_number_array(i)}) || strcmp({bond_id_new{i,2}}, c_beta{residue_number_array(i)}) || strcmp({bond_id_new{i,3}}, c_beta{residue_number_array(i)})  || strcmp({bond_id_new{i,4}}, c_beta{residue_number_array(i)}) ) ))
                    m = size(new_dihedral_parameter,1);
                    
                    for k = 1:3 %Always have 3 terms
                        new_dihedral_parameter(m+k,1) = values_tp_replace(j,k);
                        new_dihedral_parameter(m+k,2)  =  k;
                        
                        replaced = 'yes';
                        new_array_no(i) = 3;
                    end
                    
                    %Phase of cosine
                    new_dihedral_parameter(m+ 1,3) = 0;
                    new_dihedral_parameter(m+ 2,3) = 180;
                    new_dihedral_parameter(m+ 3,3) = 0;
                    
                    number_new_tp = number_new_tp + 1;
                    
                    %If replaced then move on to next dihedral value
                    j = size(values_tp_replace,1) + 1;
                    
                end
            end
            n = n + 1;
        end
        j = j + 1;
    end
    
    
    %Following part is used to differentiate between X1 and X2
    %If C alpha first then X2, this will run through all the new parameters
    %so the last replacement will be the X2 values
    if strcmp(replaced, 'yes') && (strcmp({bond_id_new{i,1}}, c_alpha{residue_number_array(i)}) || strcmp({bond_id_new{i,4}}, c_alpha{residue_number_array(i)}) )
        for l = 4:(size(names_tp_replace,1))
          
            %Forward direction
            if strcmp(OPLS_names{1}, names_tp_replace{l,1}) && strcmp(OPLS_names{2}, names_tp_replace{l,2}) && strcmp(OPLS_names{3}, names_tp_replace{l,3}) && strcmp(OPLS_names{4}, names_tp_replace{l,4})
                
                m = size(new_dihedral_parameter,1) - 3;
                
                for k = 1:3 %Always have 3 terms
                    new_dihedral_parameter(m+k,1) = values_tp_replace(l,k);
                    new_dihedral_parameter(m+k,2)  =  k;
                end
            end
            
            %Backwards direction
            if strcmp(OPLS_names{1}, names_tp_replace{l,4}) && strcmp(OPLS_names{2}, names_tp_replace{l,3}) && strcmp(OPLS_names{3}, names_tp_replace{l,2}) && strcmp(OPLS_names{4}, names_tp_replace{l,1})
                m = size(new_dihedral_parameter,1) - 3;
                
                for k = 1:3 %Always have 3 terms
                    new_dihedral_parameter(m+k,1) = values_tp_replace(l,k);
                    new_dihedral_parameter(m+k,2)  =  k;
                end
            end
            
        end
    end

     %If not replaced, use OPLS values
    if strcmp(replaced, 'no')
        n = n - 1;

        for l = 1:array_no(i)
            new_dihedral_parameter(end+1,: ) = dihedral_params(n + l,:) ;
        end
        
        n = n + 1;
    end
    
    %Moves along to next set of dihedral terms
    n = n + no;
    
    replaced = 'no';
end


%Prints to file
m = 1;
for i = 1:number_dihedrals
    for j = 1:new_array_no(i) 
        if new_dihedral_parameter(m,3) == 0
            fprintf(fid_out, '%s %s %s %s    % 1.4f   %1.0f   % 4.1f\n',  bond_id_new{i,1}, bond_id_new{i,2}, bond_id_new{i,3},  bond_id_new{i,4}, new_dihedral_parameter(m,1) , new_dihedral_parameter(m,2),  new_dihedral_parameter(m,3)  );
        else
            fprintf(fid_out, '%s %s %s %s    % 1.4f   %1.0f % 4.1f\n',  bond_id_new{i,1}, bond_id_new{i,2}, bond_id_new{i,3},  bond_id_new{i,4}, new_dihedral_parameter(m,1) , new_dihedral_parameter(m,2),  new_dihedral_parameter(m,3)  );
        end
        m = m+1;
    end
end

fclose(fid_out);
fclose(fid_log);

end
