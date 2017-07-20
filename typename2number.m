function [ numbers ] = typename2number( typenames,Atom_type )
% exchange type name to number, increase code speed and readiblity
%    
numbers=zeros(1,length(typenames));
for i_n=1:length(typenames)
    
for ix=1:length(Atom_type)
    if strcmp(typenames{i_n},Atom_type{ix})
        numbers(i_n)=ix;
        break;
    end
end

end


end

