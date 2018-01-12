function [compsOut,freqsOut] = getComponents(R1_in,R2_in,X_in,freqs_in)
R1_ind = find(freqs_in==(round(     R1_in        )),1);
H1_ind = find(freqs_in==(round(     2*R1_in      )),1);
R2_ind = find(freqs_in==(round(     R2_in        )),1);
H2_ind = find(freqs_in==(round(     2*R2_in      )),1);
E1_ind = find(freqs_in==(round(     R1_in+R2_in  )),1);
E2_ind = find(freqs_in==(round(     R1_in-R2_in  )),1);

freqsOut = [(R1_in),(2*R1_in),(R2_in),(2*R2_in),(R1_in+R2_in),(R1_in-R2_in)];
compsOut = X_in(:,[R1_ind H1_ind R2_ind H2_ind E1_ind E2_ind]);
end
