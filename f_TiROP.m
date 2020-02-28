function z = f_TiROP(x, dims, n_sur)
% Time Reversibility by Ordinal Patterns (TiROP) gives the z-value (z) of
% comparing the original divergence between a signal (x) and its reversed
% version vs, an empricial distribution of divergences that came up from the
% same proccess applied onto a set of (n_sur) surrogates of the same signal in
% all given dimensions (dims). A z-value higher than (aprox 1.9) is said to be
% irreversible by rejecting the null hypothesis.
% 
% INPUT:    x       : signal to test its time reversivility
%           dims    : dimension at with testing it (range from 3:7) if u want
%           n_sur   : number of surrogates of x
% 
% OUTPUT:   z       : z-value
% 
% JohannM
% Paris (2019)
% 
% Auxiliary files are in these toolkits
addpath /Users/johann.martinez/AA_Johann/1_Doctorate/ToolBoxes/JohannMCodes;
addpath /Users/johann.martinez/AA_Johann/1_Doctorate/ToolBoxes/Chaotic_Systems_Toolbox;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % DATA for testing this function (comment the function header and run it out)
% % clearvars; clc; close all;
% % a = load('/Users/johann.martinez/AA_Johann/3_Research/1_All_Publications/1_Papers/11_FP_test/Data_FP/CRISIS/EEGsignals_t_10b_chanFT10_T8.mat'); %   Despues
% % x = a.seizureData(1,:);     dims = 3:6;      n_sur = 50;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    div = zeros(1, numel(dims));     % original divergences
    for d = 1 : numel(dims)
        [~, ~, p] = f_Band_Pompe(x, dims(d));
        [~, ~, p_] = f_Band_Pompe(x(end:-1:1), dims(d));
        div(d) = f_divergences_JS_KL(p, p_);
    end
% figure, plot(x)
    ts_sur = IAAFT(x', n_sur);       % surrogate divergences
    div_sur = zeros(n_sur,  numel(dims));
    for sr = 1 : n_sur
        xsur = ts_sur(:,sr);
        for d = 1 : numel(dims)
            [~, ~, psur] = f_Band_Pompe(xsur, dims(d));
            [~, ~, psur_] = f_Band_Pompe(xsur(end:-1:1), dims(d));
            div_sur(sr,d) = f_divergences_JS_KL(psur, psur_);
        end
    end
    z = abs((div - mean(div_sur))./ std(div_sur));  % z-values
end