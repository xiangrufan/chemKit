function numbers = cell2num (cellofstr)
numbers=zeros(1,length(cellofstr));
for ix=1:length(cellofstr)
    numbers(ix)=str2double(cellofstr{ix});
end

end