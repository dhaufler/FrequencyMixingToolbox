function [isValid, solutionConditions] = FM_validate_triplet(r1r2_X,r1r2_Y,r1r2_Z)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

isValid = [];
solutionConditions = [];

if  checkPairwiseSolutions(r1r2_X,r1r2_Y) == 0 & ...
    checkPairwiseSolutions(r1r2_X,r1r2_Z) == 0 & ...
    checkPairwiseSolutions(r1r2_Y,r1r2_Z) == 0

    syms R1 R2 real     % assume ‘R1’ and ‘R2’ are two roots
    eqn1 = R1>0&R2>0;   % both roots > 0
    eqn2 = R1>R2;       % R1 is the larger root
    eqn3 = (r1r2_X(1)*R1+r1r2_X(2)*R2) > (r1r2_Y(1)*R1+r1r2_Y(2)*R2);   % X component greater than Y component.
    eqn4 = (r1r2_Y(1)*R1+r1r2_Y(2)*R2) > (r1r2_Z(1)*R1+r1r2_Z(2)*R2);   % Y component greater than Z component.

    SolutionDetails = solve(eqn1,eqn2,eqn3,eqn4,R1,R2,'ReturnConditions',true);  % solve for ‘x’ and ‘y’
    isValid = length(SolutionDetails.conditions);
    solutionConditions_temp = SolutionDetails.conditions;
    solutionConditions=char(subs(solutionConditions_temp,[SolutionDetails.R1,SolutionDetails.R2] ,{'R1','R2'}));
end
end

function dependent_pair = checkPairwiseSolutions(xy_I,xy_J)
syms a b integer;              % assume ‘a’, ‘b’,and ‘c’ are integers
eqn1 = xy_I(1)*a +xy_J(1)*b == 0;       % declare the equations
eqn2 = xy_I(2)*a +xy_J(2)*b == 0;       % declare the equations
eqn3 = a~=0;

SolutionDetails_Pairwise = solve(eqn1,eqn2,eqn3,[a b],'ReturnConditions',true);
dependent_pair = length(SolutionDetails_Pairwise.conditions);
end


