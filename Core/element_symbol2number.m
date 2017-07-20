 function [ number ] = element_symbol2number( symbol )
 tmp=load('symbols2number.mat');
 numbers_list=tmp.numbers;
 symbols_list=tmp.symbols;
 for iw=1:length(symbol)
 for ix =1 :length(symbols_list)
     
 if  strcmp(symbols_list{ix}, symbol{iw}) 
     number(iw)=numbers_list(ix);
 end
 
 end
 
 end
 end
 

