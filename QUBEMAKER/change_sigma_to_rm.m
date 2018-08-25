function [ rm ] = change_sigma_to_rm( sigma )
%The charm file format has a different set of parameters to the OPLS format
%this function changes the sigma from OPLS to rm charm format 
 
rm = (2 *  ( sigma^(6) ))^(1/6) /2; %division by 2 as rmin/2 in charm file 

end

