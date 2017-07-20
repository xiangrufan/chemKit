function Ebond = getEbond(bond_list,pos,Atom_info)
Ebond=0;
[NofB,~]=size(bond_list);
for ibond=1:NofB    
    Ati=bond_list(ibond,1);
    Atj=bond_list(ibond,2);
    [ r_desired ] = r_real( Ati,Atj,1,Atom_info);
    r_current=norm(pos(Ati,:)-pos(Atj,:));
    [ kij ] = getkij(Ati,Atj,Atom_info);
    Ebond=Ebond+kij*(r_current-r_desired)^2;
end


end
