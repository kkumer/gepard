# grace batch file for fitting
title "Fit to DVCS data"
subtitle ""
world 30, 0, 140, 12
#  --------  x axis ------------
xaxes scale Normal
xaxis  label "W"
#  --------  y axis ------------
yaxes scale Normal
yaxis  label " \xs\f{}"
#  --------  Lines ------------
s0 type xydy
s0 line linestyle 0
s0 line color 4
s0 symbol 1 
s0 symbol size 0.61
s0 symbol color 4   
s0 legend  "dataset5 (H1)"
#
s1 type xydy
s1 line linestyle 0
s1 line color 2
s1 symbol 2                                                                                                               
s1 symbol size 0.61
s1 symbol color 2   
s1 legend  "dataset6 (ZEUS)"
#
s2 hidden true
#
s3 hidden true
#
s4 hidden true
#
s5 hidden true
#
s6 type xy
s6 line linestyle 3
s6 line linewidth 1.5
s6 line color 4
s6 legend  "fit5"
#
s7 type xy
s7 line linestyle 1
s7 line linewidth 1.5
s7 line color 2
s7 legend  "fit6"
#
