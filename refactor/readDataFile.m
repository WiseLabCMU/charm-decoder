function [data] = readDataFile(filenameList)
%READDATAFILE Summary of this function goes here
%   Detailed explanation goes here
nFiles = size(filenameList, 1);
data = cell(1, nFiles);
% buffer = [0,0];
% read and load data from files sequentially
% note that this can be parallelized
for index = 1:nFiles
    buffer = importdata(filenameList[index]);
    data{index} = 1i*buffer(:,1) + buffer(:,2);
end