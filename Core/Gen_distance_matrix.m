function Dismatrix = Gen_distance_matrix (species,pos)
Natoms=length(species);
Dismatrix=zeros(Natoms,Natoms);

for ix =1:Natoms
    for iy =1:Natoms
        tmp=pos(ix,:)-pos(iy,:);
        Dismatrix(ix,iy)=(tmp*tmp')^0.5;
        if ix==iy
            Dismatrix(ix,iy)=999;
        end
    end
end

end