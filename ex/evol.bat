# grace batch file for evolution dependence
title ""
subtitle ""
world 9e-6, 0.001, 0.06, 1.5
legend length 12
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  label "\xx"
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
#  --------  y axis ------------
yaxes scale Logarithmic
yaxis  label "|\qsee label in legend box\Q| / |total amplitude|"
yaxis  ticklabel format general
yaxis  ticklabel prec 2
yaxis  tick major 10
yaxis  tick minor ticks 9
#  --------  Lines ------------
s0 line linestyle 1
s0 line linewidth 3.5
s0 line color 2
s0 legend  "LO predict."
#
s1 line linestyle 1
s1 line linewidth 4.5
s1 line color 2
s1 legend  "LOevol predict."
#
s2 line linestyle 2
s2 line linewidth 4.5
s2 line color 4
s2 legend  "NLOevolLO corr."
#
s3 line linestyle 3
s3 line linewidth 4.5
s3 line color 3
s3 legend  "NLOevolD corr."
#
s4 line linestyle 6
s4 line linewidth 4.5
s4 line color 13
s4 legend  "NLOevolND corr."
#
