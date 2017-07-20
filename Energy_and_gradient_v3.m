 function [ Energy, gradE ] = Energy_and_gradient_v3(species,pos,bond_list,angle_list)
%ENERGY_AND_GRADIENT 此处显示有关此函数的摘要
%   此处显示详细说明
pos=reshape(pos,length(species),3);

% parameter extraction test
% test C_R, N_R
 
load UFF_params
% [species,pos]=findgeomgjf_v3('C2H6.gjf');
for ispe=1:length(species)
    if strcmp('H',species{ispe})
        At{ispe}='H_';
    else
        At{ispe}='C_3';
    end
end
At_number=typename2number(At,Atom_type);
At_Z=cell2num(Z1(At_number));
At_r=cell2num(r1(At_number));
At_chi=cell2num(Xi(At_number));
At_theta=cell2num(theta0(At_number));
% for vdw use
At_D1=cell2num(D1(At_number));
At_x1=cell2num(x1(At_number));
% for torsion
At_Vi=cell2num(Vi(At_number));
At_Uj=cell2num(Uj(At_number));
%
Atom_info.At_Z=At_Z;
Atom_info.At_r=At_r;
Atom_info.At_chi=At_chi;
Atom_info.At_theta=At_theta;
Atom_info.At_D1=At_D1;
Atom_info.At_x1=At_x1;
Atom_info.At_Vi=At_Vi;
Atom_info.At_Uj=At_Uj;

% attaining bond list
old_geom=Geometry_v2(species,pos);
% bond_list = old_geom.get_bond_list();
if old_geom.Natoms>2
% angle_list = old_geom.get_angle_list();
[NofA,~]=size(angle_list);
end
[NofB,~]=size(bond_list);

NofAt=length(old_geom.Atoms);

%
Ebond=0;
Ebond_grad=zeros(NofAt,3);
for ibond=1:NofB
    
    Ati=bond_list(ibond,1);
    Atj=bond_list(ibond,2);
    [ r_desired ] = r_real( Ati,Atj,1,Atom_info);
    r_current=norm(pos(Ati,:)-pos(Atj,:));
    [ kij ] = getkij(Ati,Atj,Atom_info);
%     r_current-r_desired
    Ebond=Ebond+kij*(r_current-r_desired)^2;
    Ebond_grad(Ati,:)=Ebond_grad(Ati,:)+((pos(Ati,:)-pos(Atj,:))/r_current)*(2)*kij*(r_current-r_desired); 
    Ebond_grad(Atj,:)=Ebond_grad(Atj,:)+((-pos(Ati,:)+pos(Atj,:))/r_current)*(2)*kij*(r_current-r_desired); 
    
end

% Ebond


[Eangle,Eangle_grad] = getEangle(angle_list,species,pos,Atom_info);
% Extended_connectivity


Conmatrix_extend=zeros(NofAt,NofAt);
for i_atom=1:NofAt
    connected_by1bond_atoms=[old_geom.Atoms(i_atom).Neighbours.Index];
    Conmatrix_extend(i_atom,connected_by1bond_atoms)=1;
    for i_atom2 = connected_by1bond_atoms
        connected_by2bond_atoms = [old_geom.Atoms(i_atom2).Neighbours.Index];
        Conmatrix_extend(i_atom,connected_by2bond_atoms)=1;
    end
    
end


% Conmatrix_extend

[Evdw,vdw_grad]= evaluate_VDW_withGrad(pos,Conmatrix_extend,Atom_info);

Energy=Evdw+Ebond+Eangle;
gradE=vdw_grad+Eangle_grad+Ebond_grad;

% attaining angle list
gradE=reshape(gradE,[],1);

% ignore torsion, include in further version using RDKit source code &
% python
% % opt in xyz coord


 end

