function [vdwout, vdw_grad]= evaluate_VDW_withGrad(Positions,Conmatrix_extend,Atom_info)
 
At_D1=Atom_info.At_D1;
At_x1 =Atom_info.At_x1;
   [ Natoms,~]=size(Positions);
   
    vdw=0;
    vdw_grad=zeros(Natoms,3);
    for ix =1:Natoms
        for iy =1:ix
            if Conmatrix_extend(ix,iy)
%                 vdw=vdw+0;                
            else
                Di=At_D1(ix);
                Dj=At_D1(iy);
                xi=At_x1(ix);
                xj=At_x1(iy);
                r=norm(Positions(ix,:)-Positions(iy,:));
                Dij=sqrt(Di*Dj); % huge amount
                xij=sqrt(xi*xj);
                vdw  = vdw+ Dij*(-2*(xij/r)^6+(xij/r)^12);       
                forcexyz=(-12/r^2)*(Positions(ix,:)-Positions(iy,:))*Dij*(-(xij/r)^6+(xij/r)^12); % grad_E== force
                vdw_grad(ix,:)=vdw_grad(ix,:)+forcexyz; % maybe it is counted twice here ??? 
                vdw_grad(iy,:)=vdw_grad(iy,:)- forcexyz;          
             
              end
        end
    end
     vdwout=vdw;
 



end
