function [ r_realout ] = r_real( i,j,n,Atom_info)
%R_REAL 此处显示有关此函数的摘要
%   此处显示详细说明
At_r=Atom_info.At_r;
At_chi=Atom_info.At_chi;
r_en=  At_r(i)*At_r(j)* (sqrt(At_chi(i))-sqrt(At_chi(j)))^2   /  (At_chi(i)*At_r(i)+At_chi(j)*At_r(j));
r_bo=  -0.1332*(At_r(i)+At_r(j))*log(n);
 r_realout=  At_r(i)+At_r(j)+r_bo-r_en; % small error in UFF article
 
end

