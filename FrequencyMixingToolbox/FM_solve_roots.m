function [x_coef_R1,y_coef_R1,x_coef_R2,y_coef_R2] =  FM_solve_roots(ab_X,cd_Y)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
syms R1 R2 x y real;              % assume ‘a’, ‘b’,and ‘c’ are integers
eqn1 = ab_X(1)*R1 + ab_X(2)*R2 == x;       % declare the equations
eqn2 = cd_Y(1)*R1 + cd_Y(2)*R2 == y;       % declare the equations

[R1_sol,R2_sol] = solve(eqn1,eqn2,R1,R2);

syms w %note: 'children' behaves different for fractions depending on the 
%number of terms (for example 'x/2' vs 'x/2 + w'). So, I will add 'w' to
%both symbolic solutions for consistent behavior

% Get coefficients of x and y for the R1 and R2 expressions
%for R1
R1_terms = children(R1_sol+w);
%for R1 x term
R1_term_with_x = has(R1_terms,x);
if max(R1_term_with_x)==0
    x_coef_R1 = 0;
else
    R1_x_term = R1_terms(R1_term_with_x);
    x_coef_R1 = double(coeffs(R1_x_term,x));
end
%for R1 y term
R1_term_with_y = has(R1_terms,y);
if max(R1_term_with_y)==0
    y_coef_R1 = 0;
else
    R1_y_term = R1_terms(R1_term_with_y);
    y_coef_R1 = double(coeffs(R1_y_term,y));
end

%for R2
R2_terms = children(R2_sol+w);
%for R2 x term
R2_term_with_x = has(R2_terms,x);
if max(R2_term_with_x)==0
    x_coef_R2 = 0;
else
    R2_x_term = R2_terms(R2_term_with_x);
    x_coef_R2 = double(coeffs(R2_x_term,x));
end
%for R2 y term
R2_term_with_y = has(R2_terms,y);
if max(R2_term_with_y)==0
    y_coef_R2 = 0;
else
    R2_y_term = R2_terms(R2_term_with_y);
    y_coef_R2 = double(coeffs(R2_y_term,y));
end

end
