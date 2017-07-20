function [ kij ] = getkij( i,j,Atom_info)
%KIJ 此处显示有关此函数的摘要
%   此处显示详细说明
At_Z=Atom_info.At_Z;
BO=1; % 1.42 for aminoacid
kij= 664.12*(At_Z(i)*At_Z(j))/ ((r_real(i,j,BO,Atom_info))^3);
 

end

