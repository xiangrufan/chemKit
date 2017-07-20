% line search is not implemented, I uses a very small modifier in place of
% line search so convergence is horrible
% major Bug Fix. (In  V2, the e_vec update is wrong. )
clear
[species,pos_original]=findgeomgjf_v4('C2H6.gjf');
pos=reshape(pos_original,[],1);

% BFGS quasi newton 
% for loop
old_geom=Geometry_v2(species,pos_original);
bond_list = old_geom.get_bond_list();
if old_geom.Natoms>2
angle_list = old_geom.get_angle_list();
else
angle_list=[];
end




[ Energy, gradE ] =Energy_and_gradient_v3(species,pos,bond_list,angle_list);

B0=diag(ones(1,length(gradE)));
B0=B0*10000;

MaxDIISsize=8;
for ix=1:50
    if ix==1   
   pos_old=pos ;
   [Energy, gradE ] =Energy_and_gradient_v3(species,pos,bond_list,angle_list);  

   norm(gradE)
   
        dk=-B0\gradE;
       e_vec(1,:)=dk; % evec can not be linearly dependent !!, so has to do some change
       pos_vec(1,:)=pos;
       g_vec(1,:)=gradE;
        error_norm(1)=norm(dk);
        Bk=B0;  % this is CRUCIAL  ( NOt important if using correct )
gradE_old=gradE;
% x=0.003/norm(dk);
x=1;
 pos=pos+dk*x;
  else % start GDIIS  
    disp( ['opt cycle',num2str(ix)])
     
     [size_e_vec,~]=size(e_vec);
     if size_e_vec<MaxDIISsize
         DIISsize=size_e_vec;
     else
         DIISsize=MaxDIISsize;     
     end
     e_vec_tmp=e_vec;% debug
       for icol=1:DIISsize
           for irow=1:DIISsize
               DIISmatrix(icol,irow)=e_vec(size_e_vec-DIISsize+icol,:)*e_vec(size_e_vec-DIISsize+irow,:)';
           end
       end
     DIISmatrix(DIISsize+1,:)=1;
     DIISmatrix(:,DIISsize+1)=1;
     DIISmatrix(DIISsize+1,DIISsize+1)=0;  
     DIISmatrix;
     
%        DIISmatrix_tmp=DIISmatrix;
%       DIISmatrix_tmp(end,:)=DIISmatrix_tmp(end,:)+DIISmatrix_tmp(end-1,:);
%        setup.type = 'nofill';
%      [L,U] = ilu(sparse(DIISmatrix_tmp),setup);
%      lhs=(L*U)\[zeros(DIISsize,1);1];
%      c_vec= ((L*U)\DIISmatrix_tmp )\ lhs ;
% preconditionner =100*eye(DIISsize+1);
% preconditionner
%        c_vec= (  preconditionner*DIISmatrix )\ (preconditionner*[zeros(DIISsize,1);1])  ;
       c_vec= (DIISmatrix )\ ([zeros(DIISsize,1);1])  ;
%      c_vec = gmres(DIISmatrix_tmp,[zeros(DIISsize,1);1],5,1e-8,1000,L,U);
     
       c_vec=c_vec(1:end-1);
     pos_guess=pos_vec(end+1-DIISsize: end,:)'*c_vec;
     grad_guess=g_vec(end+1-DIISsize: end,:)'*c_vec;
     dk=-Bk\grad_guess;
     
     pos=pos_guess +dk;
     
     gradE_old=gradE;
     
    [ Energy, gradE ] =Energy_and_gradient_v3(species,pos,bond_list,angle_list);

   normForce= norm(gradE);
    if normForce<1e-2
        disp(['optimization finished at optcycle ', num2str(ix)])
        break
    else
        disp(['Norm foce is :', num2str(normForce)])
         disp(['ix foce is :', num2str(ix)])
    end
   
     e_vec(ix,:)=Bk\gradE; % this is wrong, need change
     pos_vec(ix,:)=pos;
     g_vec(ix,:)=gradE;
     error_norm(ix)=norm(dk);
  [size_e_vec,~]=size(e_vec);
  mean_error_global=mean(error_norm);
  
    end

end
pos_out=reshape(pos,length(species),3);

  Gen_G09_input_2016 (species,pos_out,'aaa','test','1','1')
