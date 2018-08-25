function  symettrise_and_polar_ddec( folder, N )

%The function works by first identifying the polar/symettry atoms present from the .psf file and
%then adapting the charges accordingly which are taken from the ddec.onetep
%file

%Function creates the necessary temporary files
create_temp_files(folder, N);

polar_hydrogens = {'H155','H168','H204','H240','H241','H270','H290','H301','H304','H310','H504','H513'};

molecule_names_charges = importdata(horzcat(folder, 'First_part_psf'));

%Finds the polar hydrogen atoms
list_polar = []; 
charges_original = []; 

for i=1:size(molecule_names_charges.textdata,1)
    
    polar = 0; 
    
    for j = 1:size(polar_hydrogens,2)
        if molecule_names_charges.textdata{i,6} == polar_hydrogens{j}
            polar = 1; 
            list_polar(end + 1) = i; 
        end
    end
    
    %list original charges 
    charges_original(i) = molecule_names_charges.data(i,1);
    atom_types_original{i} = molecule_names_charges.textdata{i,6}; 
    residue_id(i) = str2num(molecule_names_charges.textdata{i,3}); 
    
    polar_array(i) = polar; 
end

opls = importdata('../OPLS_Files/OPLS_sig_eps_values'); %contains OPLS values

%Find the original sig and eps value from the OPLS atom names
opls_sig = []; 
opls_eps = [];

for i = 1:size(atom_types_original,2)
       for j = 1: size(opls.textdata,1)
           if strcmp(atom_types_original{i}, opls.textdata{j})
               opls_sig(i) = opls.data(j,3);
               opls_eps(i) = opls.data(j,2);
           end
       end
end


%Finds the atoms attached to the polar hydrogen atoms

bonds = importdata(horzcat(folder, 'Bonds_psf')); 

%Sort into list 
n = 1;
r = 1;

for i=1:size(bonds,1)
    for j=1:size(bonds,2)
        if isnan(bonds(i,j)) || (bonds(i,j) > N)
           break 
        end
        if mod(n,2) == 1
            bonds_list(r,1) = bonds(i,j);
        else
            bonds_list(r,2) = bonds(i,j);
            r = r + 1;
        end
        n = n + 1;
    end
end

%Find neighbouring to polar atoms
neighbour_array = cell(1,N); 
list_neighbour = [];
for i=1:size(bonds_list,1)
    for j = 1:size(list_polar,2)
        if bonds_list(i,1) == list_polar(j)
            list_neighbour(end + 1) = bonds_list(i,2);
            neighbour_array{bonds_list(i,2)}(end + 1) = bonds_list(i,1); 
        end
        if bonds_list(i,2) == list_polar(j)
            list_neighbour(end + 1) = bonds_list(i,1);
            neighbour_array{bonds_list(i,1)}(end + 1)  = bonds_list(i,2); 
        end
    end
end

%Find groups of three hydrogens together for C, symettry will be applied to this
list_symmetry = zeros(N,3); 
for i=1:size(molecule_names_charges.textdata,1) 
    if (i + 2 < size(molecule_names_charges.textdata,1)) && strcmp(molecule_names_charges.textdata{i,6}, molecule_names_charges.textdata{(i + 1),6}) && strcmp(molecule_names_charges.textdata{i,6}, molecule_names_charges.textdata{(i + 2),6})
        list_symmetry(i, :) = [i, i+1, i + 2];
    end
end


volume_data = importdata(horzcat(folder, 'temp')); 
aim_volumes = volume_data.data(:,2) + volume_data.data(:,3); 

vfree = zeros(1, N);
bfree = zeros(1, N);
rfree = zeros(1, N);
volume_scaling = zeros(1, N); 
sig = zeros(1, N);
eps = zeros(1, N);
b_coeff = zeros(1, N);

for i = 1:size(volume_data.textdata)
    switch volume_data.textdata{i}(1)
        case 'H'
            vfree(i) = 7.6; 
            bfree(i) = 6.5; 
            rfree(i) = 1.64; 
        case 'C'
            vfree(i) = 34.4;
            bfree(i) = 46.6; 
            rfree(i) = 2.08; 
        case 'N'
            vfree(i) = 25.9;
            bfree(i) = 24.2;
            rfree(i) = 1.72; 
        case 'O'
            vfree(i) = 22.1; 
            bfree(i) = 15.6; 
            rfree(i) = 1.60; 
        case 'F'
            vfree(i) = 18.2; 
            bfree(i) = 9.5; 
            rfree(i) = 1.58; 
        case 'S'
            vfree(i) = 75.2;
            bfree(i) = 134.0;
            rfree(i) = 2.00; 
        case 'Cl'
            vfree(i) = 65.1; 
            bfree(i) = 94.6; 
            rfree(i) = 1.88; 
        otherwise
            disp('There is a problem with the atom name')
    end
    
    %Formula below can be derived from the Biomolecular Force Field...
    %paper. With conversion from A/B to sig/eps
    
    volume_scaling(i) =  aim_volumes(i) / vfree(i); 
    sig(i) = rfree(i) * (volume_scaling(i))^(1/3) * (2^(5/6));
    
    %Below is eps = Bi/2(sig^6) with 13.78 a unit conversion from Ha.Bohr^6 to kcal/mol.Ang^6
    eps(i) = ( bfree(i) * (volume_scaling(i))^2 * 13.7792544 ) / (4 * sig(i)^6 );
    
    %B values found for convenience for polar case
    b_coeff(i) = 4 * eps(i) * sig(i)^(6) ;
end

% Now apply these conditions to the charges found.

charges = importdata(horzcat(folder, 'charges.dat'));
charges = charges.data;

%If the psf and .onetep are not in the same order the atom order file can
%be used (created from zmat with reorder_ddec_values.m)

if (exist( horzcat(folder, 'atom_order'), 'file') == 2)
    atomorder = importdata('atom_order');
    new_order_eps = zeros(1, N);
    new_order_sig = zeros(1, N);
    new_order_b_coeff = zeros(1, N);
    new_order_charges = zeros(1, N);
    
    for i=1:size(atomorder,1)
        new_order_eps(atomorder(i)) = eps(i);
        new_order_sig(atomorder(i)) = sig(i);
        new_order_b_coeff(atomorder(i)) = b_coeff(i);
        new_order_charges(atomorder(i)) = charges(i);
    end
    eps = new_order_eps;
    sig = new_order_sig; 
    b_coeff = new_order_b_coeff;
    charges = new_order_charges;
end

for i = 1:size(charges,1)
    
    if isempty(neighbour_array{i}) == 0
        b_coeff_neigh = (sqrt(b_coeff(i)) + sum(sqrt(b_coeff(neighbour_array{i})))).^2;
        eps(i) = b_coeff_neigh / (4 * (sig(i))^6);
    end
    
    if i == list_symmetry(i,1) 
        %Symmetry of charges
        mean_charges = mean(charges(i:i+2));
        charges(i) = mean_charges;
        charges(i + 1) = mean_charges;
        charges(i + 2) = mean_charges;

        %Symmetry of LJ
        mean_eps = mean(eps(i:i+2));
        eps(i) = mean_eps;
        eps(i + 1) = mean_eps;
        eps(i + 2) = mean_eps;
        
        mean_sig = mean(sig(i:i+2));
        sig(i) = mean_sig;
        sig(i + 1) = mean_sig;
        sig(i + 2) = mean_sig;      
    end
    
    if polar_array(i) == 1
        %Set LJ parameters to zero
        eps(i) = 0;
        sig(i) = 0;
    end
    
end

fid = fopen(horzcat(folder,'lj_updated_eps80.dat'), 'wt');
for i = 1:N
   fprintf(fid, '%f, %f, %f \n', charges(i), sig(i), eps(i)); 
end

end