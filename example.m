clear all
clc

% UCI ionosphere data

X = load('data.txt'); % rows are observations
Y = load('y.txt');


% -----------------------------------------------------
lambda = 10^(-5); % regularization multiplier
numFeatures = 7;  % number of features

% HSCA with biased estimator
F1 = hsca0(X',Y',numFeatures,lambda); 
features1 = X*F1.V; %  F1.V - projection matrix

% HSCA with unbiased estimator
F2 = hsca1(X',Y',numFeatures,lambda);
features2 = X*F2.V; % F2.V - projection matrix


% eigenHSIC with biased estimator
F3 = eigenHsic0(X',Y',numFeatures);
features3 = X*F3.V; % F3.V - projection matrix

% eigenHsic with unbiased estimator
F4 = eigenHsic1(X',Y',numFeatures);
feautes4 = X*F4.V; % F4.V - projection matrix






