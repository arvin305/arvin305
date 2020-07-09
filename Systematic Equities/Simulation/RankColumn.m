function [ Y ] = RankColumn( X )
%RANKCOLUMN Summary of this function goes here
%   Detailed explanation goes here
[sortedValue_X , X_Ranked] = sort(X,'ascend');
[~, ~, ic] = unique(sortedValue_X);
Y(X_Ranked) = ic;
Y = Y';
end

