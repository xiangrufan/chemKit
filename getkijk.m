function [ kijk ] = getkijk( i,j,k,Rij,Rjk,TH,Atom_info )
%GETKIJK 此处显示有关此函数的摘要
%   此处显示详细说明
At_Z=Atom_info.At_Z;
Rik=sqrt(Rij^2+Rjk^2-2*Rij*Rjk*cos(TH));
kijk = 664.12*At_Z(i)*At_Z(k)*(3*Rij*Rjk*(1-cos(TH)^2)-cos(TH)*Rik^2)/Rik^5;% done,...

end

