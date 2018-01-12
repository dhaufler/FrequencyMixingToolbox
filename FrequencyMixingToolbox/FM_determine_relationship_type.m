function Relationship_Coefficients = FM_determine_relationship_type(r1r2_X,r1r2_Y,r1r2_Z)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
syms a b c integer;              % assume ‘a’, ‘b’,and ‘c’ are integers
eqn1 = r1r2_X(1)*a +r1r2_Y(1)*b +r1r2_Z(1)*c == 0;       % declare the equations
eqn2 = r1r2_X(2)*a +r1r2_Y(2)*b +r1r2_Z(2)*c == 0;       % declare the equations
eqn3 = a~=0;
eqn4 = b~=0;
eqn5 = c~=0;
[A, B, C] = solve(eqn1,eqn2,eqn3,eqn4,eqn5,[a b c]);

%for uniformity, let the first coefficient be positive
if A<0
    A=-A;
    B=-B;
    C=-C;
end

Relationship_Coefficients = [double(A) double(B) double(C)];
end

