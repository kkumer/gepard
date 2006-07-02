# grace batch file for fitting
title "Fit to DIS F\s2\N data"
subtitle ""
world 0, 0, 0.01, 1.5
#  --------  x axis ------------
xaxes scale Normal
xaxis  label "x\sBJ\N"
xaxis tick major 0.002                                                                                                   
xaxis  tick minor ticks 1
#  --------  y axis ------------
yaxes scale Normal
yaxis  label "F\s2\N"
yaxis tick major 0.5
yaxis  tick minor ticks 4
#  --------  Lines ------------
s0 type xydy
s0 line linestyle 0
s0 line color 4
s0 symbol 1                                                                                                               
s0 symbol size 0.61
s0 symbol color 4   
s0 legend  "dataset9 (H1), Q\S2\N = 8.5 GeV\S2\N"
#
s1 type xydy
s1 line linestyle 0
s1 line color 2
s1 symbol 1                                                                                                               
s1 symbol size 0.61
s1 symbol color 2   
s1 legend  "dataset0 (H1), Q\S2\N = 15 GeV\S2\N"
#
s2 hidden true
#
s3 hidden true
#
s4 type xy
s4 line linestyle 3
s4 line linewidth 1.5
s4 line color 4
s4 legend  "fit9"
#
s5 type xy
s5 line linestyle 3
s5 line linewidth 1.5
s5 line color 2
s5 legend  "fit0"
#
