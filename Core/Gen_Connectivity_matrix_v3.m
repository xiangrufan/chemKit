function Conmatrix = Gen_Connectivity_matrix_v3 (species,Dismatrix,covalent_radius)
 [uniques_symbols]=unique(species);
[ element_numbers ] = element_symbol2number( uniques_symbols );
atom_covalent_radiuses=covalent_radius(element_numbers);
for ix= 1:length(species)
    for iy = 1:length(species)
        ix_specie= species(ix);
        iy_specie= species(iy);
        ix_radius= find_element_radius (ix_specie)/100 ;
        iy_radius=find_element_radius (iy_specie)/100 ;
        Conmatrix(ix,iy)= (   Dismatrix(ix,iy)< (ix_radius+iy_radius)*1.15);
        
    end
end



    function radius = find_element_radius (specie)  
        for iq=1:length(uniques_symbols)
           if strcmp(uniques_symbols(iq),specie)
               radius=atom_covalent_radiuses(iq);
%            else
%                error('programmer error');
           end           
        end
    end
end