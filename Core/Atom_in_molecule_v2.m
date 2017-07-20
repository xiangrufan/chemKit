classdef Atom_in_molecule_v2 < Atom_in_molecule
    %ATOMINMOLECULE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties           
    end
    
    methods
        
         function myAtom_in_molecule = Atom_in_molecule_v2(symbol,Position,index,neighbours,is_unsaturated)   
             myAtom_in_molecule@Atom_in_molecule(symbol,Position,index,neighbours,is_unsaturated);
         end
         
%          function yes_or_no = Istype(thisObject, typename)
%            yes_or_no= strcmp(  thisObject.Element_symbol,typename);
%          end
         
         function neighbours_obj_out=get_Neighbours_condition(thisObject,condition)
             neighbours_obj_out_old=thisObject.Neighbours;
             iz=1;
             for ix=1:length(neighbours_obj_out_old)
                 if condition (neighbours_obj_out_old(ix))
                    neighbours_obj_out(iz)=neighbours_obj_out_old(ix);
                    iz=iz+1;
                 end                 
             end
         end
         
    end  
    methods (Static)
        
%          function yes_or_no = Istype(Atom, typename)
%            yes_or_no= strcmp(  Atom.Element_symbol,typename);
%          end
%          
    end
end

