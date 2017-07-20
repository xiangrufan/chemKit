function [Eangle,Eangle_grad] = getEangle(angle_list,species,pos,Atom_info)
old_geom=Geometry_v2(species,pos);
[NoA,~]=size(angle_list);
Eangle=0;
NofAt=length(species);
Eangle_grad=zeros(NofAt,3);
for iangle=1:NoA
% iangle=1;% tmp row

At_theta=Atom_info.At_theta;
    Ati=angle_list(iangle,1);
    Atj=angle_list(iangle,2);
    Atk=angle_list(iangle,3);
    TH_desired=(At_theta(Atj)/180)*pi;
    TH_current=old_geom.get_angle(Ati,Atj,Atk);
    Rij=old_geom.get_distance(Ati,Atj);
    Rjk=old_geom.get_distance(Atj,Atk);
    [ kijk ] = getkijk(Ati,Atj,Atk,Rij,Rjk,TH_desired,Atom_info );
    C2=1/(4*sin(TH_desired));
    C1=-4*C2*cos(TH_desired);
    C0=C2*(2* (cos(TH_desired) )^2+1);
    
    Eangle= Eangle+kijk*(  C0 + C1 * cos(TH_current) + C2 * cos(2*TH_current) );
    % copied from C++ code
 dist1=Rij;
dist2=Rjk;
    r =[ (pos(Ati,:) - pos(Atj,:))/dist1 ; (pos(Atk,:) - pos(Atj,:))/dist2 ];
    sinTheta=sin(TH_current);
cosTheta=cos(TH_current);
 dE_dTheta=-1 * kijk * (C1 * sin(TH_current) + 2.0 * C2 * sin(2*TH_current));
dCos_dS =[ (1.0 /dist1) *  (r(2,1)- cosTheta * r(1,1)),  (1.0 /dist1) *   (r(2,2) - cosTheta * r(1,2)),...
    (1.0 /dist1) *  (r(2,3) - cosTheta * r(1,3)),(1.0 /dist2) *  (r(1,1) - cosTheta * r(2,1)), ...
   (1.0 /dist2) *  (r(1,2) - cosTheta * r(2,2)),  (1.0 /dist2)* (r(1,3) - cosTheta * r(2,3))];

  Eangle_grad(Ati,1) =  Eangle_grad(Ati,1)+dE_dTheta * dCos_dS(1) / (-sinTheta);
  Eangle_grad(Ati,2) =  Eangle_grad(Ati,2) + dE_dTheta * dCos_dS(2) / (-sinTheta);
  Eangle_grad(Ati,3) =  Eangle_grad(Ati,3)+dE_dTheta * dCos_dS(3) / (-sinTheta);

  Eangle_grad(Atj,1)=  Eangle_grad(Atj,1)+dE_dTheta * (-dCos_dS(1)  - dCos_dS(4)) / (-sinTheta);
  Eangle_grad(Atj,2) = Eangle_grad(Atj,2)+ dE_dTheta * (-dCos_dS(2) - dCos_dS(5)) / (-sinTheta);
  Eangle_grad(Atj,3) = Eangle_grad(Atj,3)+dE_dTheta * (-dCos_dS(3) - dCos_dS(6)) / (-sinTheta);

  Eangle_grad(Atk,1) = Eangle_grad(Atk,1)+dE_dTheta * dCos_dS(4) / (-sinTheta);
  Eangle_grad(Atk,2) = Eangle_grad(Atk,2)+dE_dTheta * dCos_dS(5) / (-sinTheta);
  Eangle_grad(Atk,3) = Eangle_grad(Atk,3)+dE_dTheta * dCos_dS(6) / (-sinTheta);

    
    
    
end


end
