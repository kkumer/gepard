# grace batch file for plotting
# accuracy graphs
title ""
subtitle ""
world 1e-7, 1e-8, 0.9, 1
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  label "\xx"
xaxis  tick on
xaxis  tick major 10
xaxis  tick minor ticks 1
xaxis  ticklabel on
xaxis  ticklabel format power
xaxis  ticklabel prec 0
#  --------  y axis ------------
yaxes scale Logarithmic
yaxis  label "relative error"
yaxis  tick on
yaxis  tick major 10
yaxis  tick minor ticks 1
yaxis  ticklabel on
yaxis  ticklabel format power
yaxis  ticklabel prec 0
#  --------  Lines ------------
s0 type xy
s0 line linestyle 1
s0 line linewidth 1.5
s0 line color 1
s0 legend  "SPEED = 1"
#
s1 type xy
s1 line linestyle 2
s1 line linewidth 1.5
s1 line color 2
s1 legend  "SPEED = 2"
#
s2 type xy
s2 line linestyle 4
s2 line linewidth 1.5
s2 line color 4
s2 legend  "SPEED = 3"
#
s3 hidden true
#
