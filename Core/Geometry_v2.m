classdef Geometry_v2<Geometry
    %GEOMETRY_V2 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        %         Smiles;
        bond_list
    end
    
    methods
        
        function mygeom = Geometry_v2(Atoms_symb_input,Pos)
            mygeom@Geometry(Atoms_symb_input,Pos);
            load symbol2radius
            if 1==length(Atoms_symb_input)
                mygeom.Atomsymbols=Atoms_symb_input;
                Atom=Atom_in_molecule_v2(Atoms_symb_input{1},Pos,1,[],0);
                %                 mygeom.Atoms(1)
                mygeom.Atoms =Atom;
                mygeom.Natoms =1;
                mygeom.Atoms.addParent(mygeom);
            else
                Dismatrix = Gen_distance_matrix (Atoms_symb_input,Pos);
                Conmatrix= Gen_Connectivity_matrix_v3 (Atoms_symb_input, Dismatrix,covalent_radius);
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
        function bond_list = get_bond_list(mygeom)
         Conmatrix=   mygeom.Conmatrix;
         [cols,~]=size(Conmatrix);
         iz=1;
         for ix =1:cols
             for iy=1:ix
             if Conmatrix(ix,iy)==1
                 bond_list(iz,:)=[ix,iy];
                 iz=iz+1;
             end
             
             end
         end
         
         mygeom.bond_list=bond_list;
         
        
        end
        function angle_list = get_angle_list(mygeom)
            if isempty (mygeom.bond_list)
                mygeom.get_bond_list();
            end
            iz=1;
            bond_list=mygeom.bond_list;
            [cols,~]=    size(bond_list);
            for ix=1:cols
                for iy=1:cols
                    if ix~=iy
                        bond_x=bond_list(ix,:);
                        bond_y=bond_list(iy,:);
                        
                        if bond_x(1)==bond_y(2)
                            angle_list(iz,:)=[bond_x(2),bond_x(1),bond_y(1)];
                            iz=iz+1;
                        elseif  bond_x(2)==bond_y(2)
                            angle_list(iz,:)=[bond_x(1),bond_x(2),bond_y(1)];
                            iz=iz+1;
                        elseif  bond_x(1)==bond_y(1)
                            angle_list(iz,:)=[bond_x(2),bond_x(1),bond_y(2)];
                            iz=iz+1;
                        elseif  bond_x(2)==bond_y(1)
                            angle_list(iz,:)=[bond_x(1),bond_x(2),bond_y(2)];
                            iz=iz+1;
                        end
                    end                    
                end
            end
            idel=1;
            for itest1=1:(iz-1)
                 for itest2=1:itest1
                     angA= angle_list(itest1,:);
                     angB= angle_list(itest2,:);
                     if angA([3,2,1])==angB([1,2,3])
                         deletelist(idel)=itest1;
                         idel=idel+1;
                     end
                 end
            end
           angle_list( deletelist,:)=[];
        end
         function distance = get_distance(mygeom,Ati,Atj)
             distance=norm(mygeom.Atoms(Ati).Position-mygeom.Atoms(Atj).Position);
         end
            function [angle,vector1,vector2] = get_angle(mygeom,Ati,Atj,Atk)
             vector1=(mygeom.Atoms(Atj).Position-mygeom.Atoms(Ati).Position);
             vector2=(mygeom.Atoms(Atj).Position-mygeom.Atoms(Atk).Position);
             vector1=vector1/norm(vector1);
             vector2=vector2/norm(vector2);
%              (sum(vector1.*vector2)/(norm(vector1)*norm(vector2)))
             angle= acos(sum(vector1.*vector2));
             
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
                Atom=Atom_in_molecule_v2(symbol,Position,ix,neighbours,is_unsaturated)    ;
                Atoms(ix)=Atom;
            end
        end
    end
end

