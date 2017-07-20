function Conmatrix = Gen_Connectivity_matrix (species,Dismatrix)

for ix= 1:length(species)
    % THIS IS NOT SYMMETRIC
    for iy = 1:length(species)
        
    if strcmp(species(ix),'H')|| strcmp(species(iy),'H')
     Conmatrix(ix,iy)=   Dismatrix(ix,iy)<1.3;
    elseif strcmp(species(ix),'Cl')|| strcmp(species(iy),'Cl')
%         ix
%      'aaa'
%       Dismatrix(ix,iy)<1.8 
      Conmatrix(ix,iy)=   Dismatrix(ix,iy)<1.8 ;
   elseif strcmp(species(ix),'S')|| strcmp(species(iy),'S')
      Conmatrix(ix,iy)=   Dismatrix(ix,iy)<1.84;
    else
     Conmatrix(ix,iy)=   Dismatrix(ix,iy)<1.6;
    end
    
    end
end

end