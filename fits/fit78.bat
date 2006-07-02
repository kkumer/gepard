# grace batch file for fitting
title "Fit to DIS F\s2\N data"
subtitle ""
world 0, 0, 70, 1.5
#  --------  x axis ------------
xaxes scale Normal
xaxis  label "Q\S2\N [GeV\S2\N]"
xaxis tick major 20
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
s0 legend  "dataset7 (H1) x\sBJ\N = 0.002"
#
s1 type xydy
s1 line linestyle 0
s1 line color 2
s1 symbol 1                                                                                                               
s1 symbol size 0.61
s1 symbol color 2   
s1 legend  "dataset8 (H1) x\sBJ\N = 0.0005"
#
s2 type xy
s2 line linestyle 3
s2 line linewidth 1.5
s2 line color 4
s2 legend  "fit7"
#
s3 type xy
s3 line linestyle 3
s3 line linewidth 1.5
s3 line color 2
s3 legend  "fit8"
#
s4 hidden true
#
s4 hidden true
#
