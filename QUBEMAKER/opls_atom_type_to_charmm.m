function [ OPLS_names ] = opls_atom_type_to_charmm( CHARMM_name, number_atom_type )

%This script changes Charmm atom names to OPLS (want to do it this way round
%as otherwise one to many)

OPLS_names = cell(4,1);

%Go through charmm names and from the number find the corresponding OPLS
%names from the Number_to_Atom_type list.

for i=1:4
      OPLS_names{i} = number_atom_type{str2num(CHARMM_name{i}(2:end)),2};
end

end
