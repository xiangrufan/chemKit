classdef Atom_in_molecule < hgsetget
    %ATOMINMOLECULE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
           Element_symbol;
           Position;
           Index;
           Neighbours;
           Is_unsaturated;
           Parent_mol=[];
    end
    
    methods
        
         function myAtom_in_molecule = Atom_in_molecule(symbol,Position,index,neighbours,is_unsaturated)            
             myAtom_in_molecule.Is_unsaturated=is_unsaturated;
             myAtom_in_molecule.Neighbours=neighbours;
             myAtom_in_molecule.Position=Position;
             myAtom_in_molecule.Index=index;
             myAtom_in_molecule.Element_symbol=symbol;             
         end
         
         function yes_or_no = Istype(thisObject, typename)
           yes_or_no= strcmp(  thisObject.Element_symbol,typename);
         end
         
         function addParent(thisObject, parent)      
              thisObject.Parent_mol=parent;
          end
         function neighbours_obj_out=get.Neighbours(thisObject)    
% if nargin>2
%     
% end
              neighbours_obj_out=thisObject.Parent_mol.Atoms(thisObject.Neighbours);
              
         end
%          function neighbours_obj_out=get_Neighbours_condition(thisObject,condition)
%              neighbours_obj_out_old=thisObject.Neighbours;
%              iz=1;
%              for ix=1:length(neighbours_obj_out_old)
%                  if condition(neighbours_obj_out_old(ix))
%                     neighbours_obj_out(iz)=neighbours_obj_out_old(ix);
%                     iz=iz+1;
%                  end
%                  
%              end
%          end
         
    end
    methods(Static)
        
        
        
        
        function leftover= exclude (atoms,atoms2exclude)
            exlude_list=[];
            for ix=1:length(atoms)
                for iy=1:length(atoms2exclude)
%                     ix
%                     iy
%                     atoms(ix).Index
%                     atoms2exclude(iy)
% atoms(ix).Index
% atoms2exclude
%                     atoms2exclude(iy).Index
                    if atoms(ix).Index==atoms2exclude(iy).Index
                    exlude_list=[exlude_list,ix];
                    end
                end
            end
%             idx=1:length(atoms);
%             idx=idx(~idx==exlude_list)
% exlude_list
atoms(exlude_list)=[];
            leftover=atoms;% still wrong....
        end
        
        
        
        
    end
end

