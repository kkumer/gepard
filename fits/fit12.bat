# grace batch file for fitting
title "Fit to DVCS data"
subtitle ""
world 0, 0.05, 1, 100
#  --------  x axis ------------
xaxes scale Normal
xaxis  label "- t [GeV\S2\N]"
#  --------  y axis ------------
yaxes scale Logarithmic
yaxis  label "d\xs\f{}/d t"
#  --------  Lines ------------
s0 type xydy
s0 line linestyle 0
s0 line color 4
s0 symbol 1                                                                                                               
s0 symbol size 0.61
s0 symbol color 4   
s0 legend  "dataset1 (H1)"
#
s1 type xydy
s1 line linestyle 0
s1 line color 2
s1 symbol 1                                                                                                               
s1 symbol size 0.61
s1 symbol color 2   
s1 legend  "dataset2 (H1)"
#
s2 type xy
s2 line linestyle 3
s2 line linewidth 1.5
s2 line color 4
s2 legend  "fit1"
#
s3 type xy
s3 line linestyle 1
s3 line linewidth 1.5
s3 line color 2
s3 legend  "fit2"
#
s4 hidden true
#
s5 hidden true
#
s6 hidden true
#
s7 hidden true
#
