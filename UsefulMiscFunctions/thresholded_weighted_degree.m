function [wdeg_vector] = thresholded_weighted_degree(X, threshold)
% This function takes in a diagonally symmetrical non-negative matrix and
% computes thresholded weighted node degree for each node.
%
% By: Diego G. Davila 
%     Proekt Lab 
%     University of Pennsylvania School of Medicine
%     10/29/2021
% 
% INPUTS: 
%     1. X: A diagonally symmetrical non-negative n x n matrix. 
% 
% OUTPUTS:
%     1. wdeg_vector: A vector of size n, where each element is the sum of elements along the rows of X that exceed a threshold value(thresholded weighted degree).
%
% ----------------------------------------------------------------------------------------------------------------------------------
% First, make sure the diagonals are zero (to make summing procedure easier)
X = X - diag(diag(X));
% set value sbelow threshold to zero
X(X < threshold) = 0;
% now, sum along the rows
wdeg_vector = sum(X, 2);