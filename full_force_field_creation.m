% Main Function for creating NAMD input files
% for protein specific force field. Run this script for file creation

% Openmm file input specified by changing variable to Y
openmm = 'N';

%Name of folder with protein files
folder_name = 'Example_1P7E';

folder = horzcat(  '/', folder_name, '/');
inputfolder  = horzcat('../Input_File', folder , '');

%Name of backbone terms used for ala/gly backbone dipeptides
tp_name= 'backbone_tp_0.50';

%Number of atoms
text = fileread(horzcat(inputfolder,'ddec.onetep'));
numb_ind = strfind(text, 'Totals');
N = str2double(text((numb_ind+9):(numb_ind+17)));

symettrise_and_polar_ddec( inputfolder, N );

lj_params = importdata(horzcat(inputfolder,'lj_updated_eps80.dat'));

%Make output folders 
if exist( horzcat('../Output_File/',folder_name), 'dir') == 0 
    mkdir(horzcat('../Output_File/',folder_name))
end

if exist( horzcat('../Final_File/',folder_name), 'dir') == 0 
    mkdir(horzcat('../Final_File/',folder_name))
end

%Find the masses of the atoms
get_mass( folder, N );

%Get bond/angles/dihedrals/impropers and new ionised file  
[sum_charge_before, sum_charge_after] = process_ionised_psf( inputfolder, N );
copyfile(horzcat(inputfolder,'new_ionized.psf'), horzcat('../Output_File',folder));
copyfile(horzcat(inputfolder,'new_ionized.psf'), horzcat('../Final_File/',folder));

%Angles and bonds
script_getbondedparams( folder, N )
script_angleparams( folder, N );

%LJ parameters
script_LJ( folder, openmm );
 
%Improper
script_angleimproper( folder, N );

%Dihedrals
script_angledihedral( folder, N, tp_name );

%Puts all the force field components together in a file in the Final_File folder
final_force_field( folder, tp_name, openmm );
