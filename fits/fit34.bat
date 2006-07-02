# grace batch file for fitting
title "Fit to DVCS data"
subtitle ""
world 0, 0.02, 90, 25
#  --------  x axis ------------
xaxes scale Normal
xaxis  label "-q\s1\N\S2\N"
#  --------  y axis ------------
yaxes scale Logarithmic
yaxis  label " \xs\f{}"
#  --------  Lines ------------
s0 type xydy
s0 line linestyle 0
s0 line color 4
s0 symbol 1
s0 symbol size 0.61
s0 symbol color 4   
s0 legend  "dataset3 (H1)"
#
s1 type xydy
s1 line linestyle 0
s1 line color 2
s1 symbol 2
s1 symbol size 0.61
s1 symbol color 2   
s1 legend  "dataset4 (ZEUS)"
#
s2 hidden true
#
s3 hidden true
#
s4 type xy
s4 line linestyle 3
s4 line linewidth 1.5
s4 line color 4
s4 legend  "fit3"
#
s5 type xy
s5 line linestyle 1
s5 line linewidth 1.5
s5 line color 2
s5 legend  "fit4"
#
s6 hidden true
#
s7 hidden true
#
