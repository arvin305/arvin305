function [Pos] = ReturnColumnNumber(T,StringPattern)
%RETURNCOLUMNNUMBER Summary of this function goes here
%   Detailed explanation goes here

[m,n] = size(T);
Pos = [];
for i = 1:n
    bool = strfind(T.Properties.VariableNames{i},StringPattern)  ;  
    if bool > 0
        Pos = [Pos,i];
    end
end

