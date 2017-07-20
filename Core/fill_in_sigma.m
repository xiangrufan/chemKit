function Electron_needed_list= fill_in_sigma(Atoms,numberofbond_list)
Valence_Electron_dist=Gen_Valence_Electron_list(Atoms);
 [ elem_numbers ] = element_symbol2number( Atoms );
 
sigmabond_Valence_Electron_dist=Valence_Electron_dist+numberofbond_list;
for ix =1:length(Atoms)
   elem_number= elem_numbers(ix);
 if elem_number<2
    Electron_needed_list(ix)= 2-sigmabond_Valence_Electron_dist(ix);% the electron needed for 
    %full shell  minus electron it have, including share
 elseif elem_number<18
     Electron_needed_list(ix)=8-sigmabond_Valence_Electron_dist(ix);
 elseif elem_number<54
     Electron_needed_list(ix)=18-sigmabond_Valence_Electron_dist(ix);
 else
%       disp('unsupported_element, ANALYSIS module will be wrong');
     Electron_needed_list(ix)=0;
 end    
end 


end
