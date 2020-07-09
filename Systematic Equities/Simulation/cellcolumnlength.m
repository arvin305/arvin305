function [length] =  cellcolumnlength(A)
    length = min(find(strcmp(A,''))) - 1;
end