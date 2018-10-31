function  script_LJ( folder, openmm )

%LJ parameters

%This script gets the C6 and C12 terms used in the force field from the
%DDEC values and assigns them there atom type. This is where 0.5,0.5 or 1,1 LJ scaling is decided
% onefourscaling decides is the 1-4 LJ interaction should be scaled by 0.5
% as happens in the OPLS force field 

%Input
inputfolder  = horzcat('../Input_File', folder , '/');

%Output
fid = fopen(horzcat('../Output_File',folder,'lj_ddec'), 'w');
new_psf = importdata(horzcat(inputfolder,'new_AA_psf'));

lj_params = importdata(horzcat(inputfolder,'lj_updated_eps80.dat'));

sigma = lj_params(:, 2);
epsilon = lj_params(:, 3);

if ( strcmpi((openmm), 'Y') == 0 ) 
    epsilon = - epsilon; %Charm convention
end

rmin_over_two = zeros(size(lj_params,1),1);

for i = 1:size(lj_params,1)
    %Change format to that used in charmm file format
    rmin_over_two(i) = change_sigma_to_rm(sigma(i));
end

for i = 1:size(lj_params,1)
    %To use 0.5,0.5 scaling change the final epsilon to be
    %divided by two
    fprintf(fid, '%s  %1.2f  %1.6f   %1.6f  %1.2f  %1.6f   %1.6f \n', new_psf.textdata{i, 6}, 0.00, epsilon(i), rmin_over_two(i), 0.00, (epsilon(i))/2, (rmin_over_two(i)) );
end

end