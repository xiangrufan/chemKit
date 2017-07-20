[species,pos_original]=findgeomgjf_v4('C2H6.gjf');
pos=reshape(pos_original,[],1)


% for loop
old_geom=Geometry_v2(species,pos_original);
bond_list = old_geom.get_bond_list();
if old_geom.Natoms>2
angle_list = old_geom.get_angle_list();
else
angle_list=[];
end



for ix=1:100
    pos_old=pos;
    gradE_old=gradE;
  [Energy, gradE ] =Energy_and_gradient_v3(species,pos,bond_list,angle_list);      
%   norm(gradE) 
  if norm(gradE) >0.03 && norm(gradE_old) >0.03% no use... gradE always >0.03 
%        modifier=0.005/norm(gradE) % unstable..... zig-zagging
       % modifier has to be extremely small
    modifier=1/10000; %  stable..... 
  else
        modifier=1/10000;
%         modifier=1;
  end
  
  pos=pos-modifier*gradE   ; 
  
   norm(gradE)
end
pos
pos_out=reshape(pos,length(species),3);
%%
  Gen_G09_input_2016 (species,pos_out,'aaa','test','1','1')