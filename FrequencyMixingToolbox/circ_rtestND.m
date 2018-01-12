function [pval z] = circ_rtestND(alpha, dim)
%
% MODIFIED from 'circ_rtest' by Darrell Haufler, 8/15/2017 for use with
% Frequency Mixing Toolbox functions. Original documentation for circ_rtest
% below:
% [pval, z] = circ_rtest(alpha,w)
%   Computes Rayleigh test for non-uniformity of circular data.
%   H0: the population is uniformly distributed around the circle
%   HA: the populatoin is not distributed uniformly around the circle
%   Assumption: the distribution has maximally one mode and the data is 
%   sampled from a von Mises distribution!
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%
%   Output:
%     pval  p-value of Rayleigh's test
%     z     value of the z-statistic
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


n = length(alpha);
r =  circ_r(alpha,[],[],dim);

% compute Rayleigh's R (equ. 27.1)
R = n.*r;

% compute Rayleigh's z (equ. 27.2)
z = R.^2 / n;

% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4*n+4*(n^2-R.^2))-(1+2*n));








