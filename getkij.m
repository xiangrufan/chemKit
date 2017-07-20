function [ kij ] = getkij( i,j,Atom_info)
%KIJ �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
At_Z=Atom_info.At_Z;
BO=1; % 1.42 for aminoacid
kij= 664.12*(At_Z(i)*At_Z(j))/ ((r_real(i,j,BO,Atom_info))^3);
 

end

