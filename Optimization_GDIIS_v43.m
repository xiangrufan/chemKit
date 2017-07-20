% line search is not implemented, I uses a very small modifier in place of
% line search so convergence is horrible
% V3 major Bug Fix. (In  V2, the e_vec update is wrong. )
% V4 using BFGS hessian update,  
% V41 improve BFGS hessian update,  implemented Trustradius 
% V42 exactt hessian
clear
[species,pos_original]=findgeomgjf_v4('C2H6.gjf');
pos=reshape(pos_original,[],1);


old_geom=Geometry_v2(species,pos_original);
bond_list = old_geom.get_bond_list();
if old_geom.Natoms>2
    angle_list = old_geom.get_angle_list();
else
    angle_list=[];
end



[ Energy, gradE ] =Energy_and_gradient_v3(species,pos,bond_list,angle_list);

B0=diag(ones(1,length(gradE)));
B0=B0*1000;
%      B0 = Gen_Hessian (@(posX) Energy_and_gradient_v3(species,posX,bond_list,angle_list),pos);

MaxDIISsize=12;
trustradius=0.03;
for ix=1:31
    if ix==1   
 
   [Energy, gradE ] =Energy_and_gradient_v3(species,pos,bond_list,angle_list);  
   if norm(gradE)>10000;
       gradE=100*gradE/norm(gradE);
   end
       
       
   norm(gradE)
   % trust radius   
        dk=-B0\gradE;
        if norm(dk)>trustradius
         modifier=  trustradius/norm(dk)
         dk=dk*modifier;       
        end
%      eig(B0)
       e_vec(:,1)=-B0\gradE ;% evec can not be linearly dependent !!, so has to do some change this is Bug, done
       pos_vec(:,1)=pos;
       g_vec(:,1)=gradE;
       Energy_vec(1)=Energy;
        error_norm(1)=norm(gradE);
        Bk=B0;   

 

 pos=pos+dk;
  else % start GDIIS  
      
    disp( ['opt cycle ',num2str(ix)])
     
     [~,size_e_vec]=size(e_vec);
     if size_e_vec<MaxDIISsize
         DIISsize=size_e_vec;
     else
         DIISsize=MaxDIISsize;     
     end
     e_vec_tmp=e_vec;% debug
       for icol=1:DIISsize
           for irow=1:DIISsize
               DIISmatrix(icol,irow)=e_vec(:,size_e_vec-DIISsize+irow)'*e_vec(:,size_e_vec-DIISsize+icol);
           end
       end
     DIISmatrix(DIISsize+1,:)=1;
     DIISmatrix(:,DIISsize+1)=1;
     DIISmatrix(DIISsize+1,DIISsize+1)=0;  
     DIISmatrix;
     
 
       c_vec= (DIISmatrix )\ ([zeros(DIISsize,1);1])  ;
 
     
       c_vec=c_vec(1:end-1);
     pos_guess=pos_vec(:,end+1-DIISsize: end)*c_vec;
     grad_guess=g_vec(:,end+1-DIISsize: end)*c_vec;
       % this is not updated, strange
       dk_old=-(Bk+100*eye(length(Bk)))\grad_guess;
     dk=dk_old;
    
  Energy_guess=Energy_vec(end+1-DIISsize: end)*c_vec;
     
    if norm(dk)>trustradius
        dk=dk*trustradius/norm(dk); % aleady good              
    end
      
      
    
     pos_old=pos  ;
     pos=pos_guess +dk;     
     gradE_old=gradE;     
    [ Energy, gradE ] =Energy_and_gradient_v3(species,pos,bond_list,angle_list);
%     if norm(gradE)>10000;
%        gradE=100*gradE/norm(gradE);
%    end
    
%     Energy_pre=dk*grad_guess
% Energy_guess 
% Energy_guess=Energy_and_gradient_v3(species,pos_guess,bond_list,angle_list)
% use energy_old not energy_guess, avoid error. 
      ratio= (Energy-Energy_vec(end))/(g_vec(:,end)'*dk+0.5*dk'*Bk*dk);

    
 if ratio<0 || (Energy-Energy_vec(end) )>0% need seprate function , reject current step, recalc energy, redo 
     disp('Wrong direction of search, Using BFGS instead')
        trustradius=trustradius*0.25;    
        
%         eig(Bk_new)
        ratio_before_adjust=ratio
       Is_energy_increase=  (Energy-Energy_vec(end) )>0
% Bk = Gen_Hessian (@(posX) Energy_and_gradient_v3(species,posX,bond_list,angle_list),pos);
Bk = force_symmetry(Bk);
% trustradius=0.03;% not need to restrict too much
            dk=-(Bk+100*eye(length(Bk)))\g_vec(:,end);
%         dk=(grad_guess);
%        dk=-Bk\grad_guess;
if dk>trustradius
        dk=dk*trustradius/norm(dk); % aleady good              
end

     pos=pos_vec(:,ix-1) +dk; 
     [ Energy, gradE ] =Energy_and_gradient_v3(species,pos,bond_list,angle_list);
      ratio= (Energy-Energy_vec(end))/(g_vec(:,end)'*dk+0.5*dk'*Bk*dk)
      Is_energy_increae_after_adjust=  (Energy-Energy_vec(end) )>0
      
      
%     ratio= (Energy-Energy_guess)/(grad_guess'*dk+0.5*dk'*Bk*dk);
 elseif ratio<0.25
      trustradius=trustradius*0.25;    
 elseif  ratio>0.75
     trustradius=trustradius*2;
 end
    
 if trustradius<1e-8
     trustradius=1e-8;
 elseif trustradius>0.03
      trustradius=0.03;
 end
 trustradius
 % Trying true Hessian as input. This way we can test mettle of BFGS
%        yk=gradE-grad_guess; %  BFGS Update makes solution worse         
%        sk=pos-pos_guess; % check this.
if ix>2
    
%     yk=g_vec(:,ix-1)-g_vec(:,ix-2);
%     sk=pos_vec(:,ix-1)-pos_vec(:,ix-2);
     yk=gradE-g_vec(:,ix-1);
    sk=pos-pos_vec(:,ix-1);
    tk=yk-Bk*sk;
    if all(norm(sk)<0.001)
       disp( 'Not updated due to linear dependence');
        Bk_new=Bk;
    else
%            Bk_old
       Bk_new=Bk_old+(yk*yk')/(yk'*sk)-(Bk_old*(sk*sk')*Bk_old)/ (sk'*Bk_old*sk);% why this is always slower???
%    Bk_new=Bk+ (yk*yk')/(yk'*sk)-(Bk*sk*sk'*Bk')/ (sk'*Bk*sk);  % seems wrong this update
%                Bk_new=Bk+(tk*sk'+sk*tk'-sk*((tk'*sk)/(sk'*sk))*sk')*inv(sk'*sk);
        %     Bk_new=Bk;
    end
else
      Bk_new=Bk;
end
     
    
    
   normForce= norm(gradE)
    if normForce<1e-2
        disp(['optimization finished at optcycle ', num2str(ix)])
        break
    else
        disp(['Norm foce is :', num2str(normForce)])
         disp(['ix foce is :', num2str(ix)])
    end
   
     e_vec(:,ix)=-Bk\gradE; 
     pos_vec(:,ix)=pos;
     g_vec(:,ix)=gradE;
     Energy_vec(ix)=Energy;
     Bk_old=Bk;
    error_norm(ix)=norm(gradE);
  [size_e_vec,~]=size(e_vec);

for iaa=1:length(Bk) % remove noise, Bk should be symmetric
    for ibb=1:length(Bk)
        Bk_new(iaa,ibb)=Bk_new(ibb,iaa);
    end
end

  [V,D] = eig(Bk_new);
  [~,loe]=size(D);
  Bk=Bk_new; % slower..... WHY ...
 eigs1= eig(Bk);  
    end
end

pos_out=reshape(pos,length(species),3);

Gen_G09_input_2016 (species,pos_out,'aaa','test','1','1')
