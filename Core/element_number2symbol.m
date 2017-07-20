 function [ symbol ] = element_number2symbol( number )
 tmp=load('symbols2number.mat');
 numbers_list=tmp.numbers;
 symbols_list=tmp.symbols;
 for iw=1:length(number)
%  for ix =1 :length(symbols_list)
     
symbol{iw}=symbols_list{number(iw)};
 
%  end
 
 end
 end
 

