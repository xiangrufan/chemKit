classdef Geometry< hgsetget
    %GEOMETRY Summary of this class goes here
    %   Detailed explanation goes here
    

    properties
        Atoms;% THIS IS ATOM OBJECTS
        Atomsymbols;
        Pos;% specially implemented! 
        Natoms;     
        Dismatrix;
        Conmatrix;     
    end
    
methods
        function mygeom = Geometry(Atoms_symb_input,Pos)       
            if 1==length(Atoms_symb_input)              
                mygeom.Atomsymbols=Atoms_symb_input;
                Atom=Atom_in_molecule(Atoms_symb_input{1},Pos,1,[],0);
%                 mygeom.Atoms(1)
                mygeom.Atoms =Atom;   
                mygeom.Natoms =1;   
                mygeom.Atoms.addParent(mygeom);
            else
                    Dismatrix = Gen_distance_matrix (Atoms_symb_input,Pos);
                    Conmatrix= Gen_Connectivity_matrix (Atoms_symb_input, Dismatrix);             
                    mygeom.Conmatrix=Conmatrix;
                    mygeom.Dismatrix=Dismatrix;
                    mygeom.Atomsymbols=Atoms_symb_input;
        %             mygeom.Pos=Pos;
                    mygeom.Natoms= length(Atoms_symb_input);
                    mygeom.Atoms=mygeom.assign_atoms(Atoms_symb_input,Pos,Conmatrix);
                    for ix=1: length(mygeom.Atoms)
                        mygeom.Atoms(ix).addParent(mygeom);
                    end
            end
        end    
        function already_selected = select_connected_atoms(mygeom,node2be_searched)                             
                    conmatrix=mygeom.Conmatrix;
                    fast_select_neighbour= @(current_selected) find(conmatrix(current_selected,:)==1); 
                    already_selected=[];
                       iy=1;
                    while ~isempty(node2be_searched)
                        iy=iy+1;
                     if iy>500
                        disp( 'Algorithm complexity limit reached, aborted to save time');
                         break     
                     end
                        already_selected=  union (node2be_searched,already_selected);  
                        next_selected=[];
                        for ix=1:length(node2be_searched)
                        tmp=fast_select_neighbour(node2be_searched(ix)) ;  % this function ONLY ACCEPT single input, ARRAY input gives WRONG output
                        next_selected=union(tmp,next_selected);
                        end
                       node2be_searched= setdiff(next_selected,already_selected);

                    end
        end
        function output_Pos=get.Pos(mygeom)    
            output_Pos=zeros(mygeom.Natoms,3);
          for ix=1:mygeom.Natoms
         output_Pos(ix,:)=   mygeom.Atoms(ix).Position;
          end
        
        end       
        function Atom = get_atom(mygeom,atom_NO)
            if isnumeric(atom_NO)
                Atom=mygeom.Atoms(atom_NO);
            elseif ischar(atom_NO)
                for ix=1:length(mygeom.Atomsymbols)
                    if atom_NO==mygeom.Atomsymbols{ix}
                        Atom=mygeom.Atoms(ix);
                        break;
                    end
                end
            end
        end
        
         function atom_pos = get_atom_pos(mygeom,atom_NO)
             % output position in [x,y,z] format
         Atom=mygeom.get_atom(atom_NO);
             atom_pos=Atom.Position;
         end               
         function gen_moleculeplot(mygeom)
        fid4= fopen(['temp/tmp.smi'],'w');
          fprintf(fid4, '%s',mygeom.Smiles) ;
         fclose(fid4);
         dos(['obabel temp/tmp.smi -O  temp/',mygeom.Unique_name,'.png']);
         end
         function multiplicity=get_multiplicity(mygeom,charge)
             total_electrons =   mygeom.get_total_electrons(charge);
             if 1==rem(total_electrons,2)
                 multiplicity=2;
             else
                 multiplicity=1;
             end
         end
         function [species, pos]= output_geom(this)
             species=    this.Atomsymbols;
             pos      =    this.Pos;
         end
         function [num_of_X_atom]= get_num_of_X_atom(this,atomtype)
             atomlist=   this.Atomsymbols;
             num_of_X_atom=0;
             for ix=1:length(atomlist)
                 if strcmp(atomlist{ix},atomtype)
                     num_of_X_atom=num_of_X_atom+1;
                 end
             end
         end
         function [total_electrons]= get_total_electrons(this,charge)
             %            total_electrons=  this.Atomsymbols;
             [ numbers ] = element_symbol2number( this.Atomsymbols );
             total_electrons=sum(numbers)+charge;
         end
end
     
methods(Static)       
     function Atoms=assign_atoms(Atomsymbols,Pos,Conmatrix)
      
         numberofbond_list =   sum(Conmatrix );
         Electron_needed_list= fill_in_sigma(Atomsymbols,numberofbond_list);
            for ix=1:length(Atomsymbols)
             neighbours=find(Conmatrix(ix,:)==1);
             Position=Pos(ix,:);
             symbol=Atomsymbols{ix}; 
             is_unsaturated= Electron_needed_list(ix)>0;
             Atom=Atom_in_molecule(symbol,Position,ix,neighbours,is_unsaturated)    ;
             Atoms(ix)=Atom;
            end
     end
     function combined_geometry = combine_geometries (geometry1,geometry2)
        combined_geometry=Geometry( [geometry1.Atomsymbols,geometry2.Atomsymbols], [geometry1.Pos;geometry2.Pos]); 
     end
end  
 
    
end

