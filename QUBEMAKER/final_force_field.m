%Takes the new bond, angle, ddec, dihedral terms and adds them to the charm
%force field
function  final_force_field( folder, tp_name, openmm )

%Inputs
input_file_folder = horzcat('../Output_File',folder);

if strcmpi(openmm, 'Y') 
    fid = fopen(horzcat('../Final_File/inp_file_intro_openmm.inp'), 'r'); %Boss types Zmat is fine
else
    fid = fopen(horzcat('../Final_File/inp_file_intro.inp'), 'r'); %Boss types Zmat is fine  
end


%Outputs
outputfilefolder  = horzcat('../Final_File', folder);
outputfile = horzcat('./QUBE_FF', '.inp');
fid_output = fopen(outputfile, 'w'); %Boss types Zmat is fine

tline = fgets(fid);

setenv('output',outputfilefolder);
setenv('input',input_file_folder);
setenv('tp_name',tp_name);

%Scrolls through Michael's OPLS to Charmm file and adds new
%bonds,angles,dihedrals,improper and LJ in the appropriate place. Remember
%that charges are stored in the ionised.psf file 
while ischar(tline) 
        
    if strcmp(strtrim(tline), '!END') && strcmpi(openmm, 'Y')
        fclose(fid_output);
        ! cat "${input}masses" >> "QUBE_FF.inp"
        fid_output = fopen(outputfile, 'a'); %Boss types Zmat is fine
    end
    
    if strcmp(strtrim(tline), 'ANGLE')
        fclose(fid_output);
        ! cat "${input}bonds_ddec" >> "QUBE_FF.inp"
        fid_output = fopen(outputfile, 'a'); %Boss types Zmat is fine
    end
    
    if strcmp(strtrim(tline), 'DIHEDRAL')
        fclose(fid_output);
        ! cat "${input}angles_ddec" >> "QUBE_FF.inp"
        fid_output = fopen(outputfile, 'a'); %Boss types Zmat is fine
    end
    
    if strcmp(strtrim(tline), 'IMPROPER')
        fclose(fid_output);
        ! cat "${input}dihedral_ddec_${tp_name}" >> "QUBE_FF.inp"
        fid_output = fopen(outputfile, 'a'); %Boss types Zmat is fine
    end
    
    if strcmp(strtrim(tline), 'NONBONDED nbxmod 5 atom cdiel switch vatom vdistance vswitch -')
        fclose(fid_output);
        ! cat "${input}improper_ddec" >> "QUBE_FF.inp"
        fid_output = fopen(outputfile, 'a'); %Boss types Zmat is fine
    end
    
    if strcmp(strtrim(tline), 'HBOND CUTHB 0.5')
        fclose(fid_output);
        ! cat "${input}lj_ddec" >> "QUBE_FF.inp"
        fid_output = fopen(outputfile, 'a'); %Boss types Zmat is fine
    end

    fprintf(fid_output, '%s', tline);
    tline = fgets(fid);
    
end

fclose(fid);

%Move file to correct folder
if  strcmpi(openmm, 'Y')
    movefile('./QUBE_FF.inp', './QUBE_FF.par');
    movefile(horzcat('./QUBE_FF.par'), outputfilefolder);
else
    movefile(horzcat('./QUBE_FF.inp'), outputfilefolder);
end

end