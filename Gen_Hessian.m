function Hessian = Gen_Hessian (fun,variable_list)
delta=0.005;
old_list=variable_list;
for ix=1:length(variable_list)    
    new_list=old_list;
    new_list(ix)=new_list(ix)+delta/2;
   [~,grad1]= fun(new_list);
    new_list=old_list;
    new_list(ix)=new_list(ix)-delta/2;
   [~,grad2]= fun(new_list);
   
    Hessian(:,ix)=(grad1-grad2)/delta;
    
end

end