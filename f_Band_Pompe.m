function [Fp, Op, p_ord, permutations] = f_Band_Pompe(ts_1, D)
% % Function for getting the Band & Pompe parameters
% INPUTS:
% TS_1: time series
% D: embadding dimension
% OUTPUTS:
% FP: Fobidden patterns of a time series
% Op: Available patterns for the time series
% p_ord: probability of ordinal patterns
% permutacions: all posibles ordinal patterns
%
% JohannM
% Lille (2017)
%
% % % % This part is for testing purposes
% % % % clearvars; clc; close all;
% % % % load('/Volumes/JohannM/AA_Johann/3_Research/11_DataMEG/allmatrices/100307_MEG_3-Restin_rmegpreproc.mat');
% % % % ts_1 = data.trial{1}(23,:);
% % % % D = 6;                                              %Dimension of Ordinal Patterns (OrdPat)
% % % % % plot(ts_1)

if factorial(D) >= numel(ts_1)           % Condition for NOT working with Band and Pompe
    error('Take care, M <= D!  The condition must be: M >= (D+1)!, or as least M>=D! but careful with this last one');
else                                    % D! << M, always
    patterns = (numel(ts_1) - D) + 1;                             %Number of possible partitions (D-based) of my time serie with 1 sample step
    ord_Patterns = zeros(patterns, D);                  %Ordinal Patterns whose appear in my time series (D- permutations based)
    permutations = perms(0 : D - 1);                    %All possible OrdPats based dimmension D
    c = 0;  %auxiliar index for final sample
    for op = 1 : patterns                               %along all posible partitions (D-dependency)
        [~, id] = sort(ts_1(op : D + c));               %, id = previous indexes of new sortered values
        c = c + 1;                                      %increasing my auxiliar index
        
        %%%%% This part is for achieve the BP method. Associated D-ordinal pattern (0, 1, ..., D-1)
        ordinalP = zeros(1, D);                         % One specific ordinal pattern based one specific parition
        ordp = 0;       %first element in OrdPat is always zero
        for idx = 1 : D                             %along all elements in one partition of dimmension D
            ordinalP(id(idx)) = ordp;               %creating my new OrdPat according to values in each aprtition
            ordp = ordp + 1;                        %increasing ordinal pattern element until (D - 1)
        end
        ord_Patterns(op, :) = ordinalP;             %storing ONE at time OrdPat whose appear in time series
    end
    
    %%%%%% This part is for creating the probability distribution P of appearing a specific OrdPat in time serie
    [~,Op] = ismember( ord_Patterns , permutations , 'rows');                       %Ordinal patterns indexes per time serie
    p_ord = histcounts(Op, 1:factorial(D)+1, 'Normalization', 'probability');   %probability distribution P of different ordinal patterns
    Fp = find(p_ord==0);                                                            %Forbiden patterns indexes per time serie
end