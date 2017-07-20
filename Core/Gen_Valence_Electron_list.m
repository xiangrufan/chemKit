function Valence_Electron_dist = Gen_Valence_Electron_list (Atoms)

 [ numbers ] = element_symbol2number( Atoms );
for ix =1:length(Atoms)
   elem_number= numbers(ix);
 if elem_number<2
     Valence_Electron_dist(ix)=1;
 elseif elem_number<10
       Valence_Electron_dist(ix)=elem_number-2;
 elseif elem_number<18
      Valence_Electron_dist(ix)= elem_number-10;
       elseif elem_number<36
      Valence_Electron_dist(ix)= elem_number-18;
        elseif elem_number<54
      Valence_Electron_dist(ix)= elem_number-36;
      % lanthnide is special not included
 end
    
end

end