# grace batch file for evolut
page size 500 350
page background fill off
arrange(1, 2, .1, .2, .1)
with string
string on
string char size 0.780000
string 0.44, 0.80
string def "q=qS"
with string
string on
string char size 0.780000
string 0.96, 0.73
string def "q=qG"
##############################################
with g0
##############################################
world 1e-4, 0, 0.01, 4.0
view 0.15, 0.45, 0.6, 0.85
legend loctype view
legend 0.20, 0.62
legend length 5
legend char size 0.64
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
xaxis  label "iks"
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
#  --------  y axis ------------
yaxes scale Normal
yaxis  label "slopeB"
yaxis  label place spec
yaxis  label place 0.07, 0.08
yaxis  ticklabel format decimal
yaxis  ticklabel prec 0
yaxis  tick major 1.0
yaxis  tick minor ticks 4
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
#  --------  Lines ------------
s0 line pattern 0
s0 fill type 1
s0 fill rule 0
s0 fill color 7
#
s1 hidden true
s1 line linestyle 2
s1 line color 1
s1 line linewidth 2.5
s1 legend  "      LO"
#
s2 line linestyle 5
s2 line color 15
s2 line linewidth 2.5
s2 legend  "   NLO, \oMS\O"
#
s3 line linestyle 4
s3 line color 4
s3 line linewidth 2.5
s3 legend  "   NLO, \oCS\O"
#
s4 line linestyle 1
s4 line color 2
s4 line linewidth 2.5
s4 legend  "NNLO, \oCS\O"
#
##############################################
with g1
##############################################
world 1e-4, 0, 0.01, 4.0
view 0.68, 0.45, 1.13, 0.85
legend loctype view
legend 0.94, 0.87
legend length 5
legend char size 0.64
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
xaxis  label "iks"
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
#  --------  y axis ------------
yaxes scale Normal
#yaxis  label "XPDF"
#yaxis  label place spec
#yaxis  label place 0.02, 0.11
yaxis  ticklabel format decimal
yaxis  ticklabel prec 0
yaxis  ticklabel prec 0
yaxis  tick major 1
yaxis  tick minor ticks 4
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
#  --------  Lines ------------
s0 line pattern 0
s0 fill type 1
s0 fill rule 0
s0 fill color 7
#
s1 hidden true
s1 line linestyle 2
s1 line color 1
s1 line linewidth 2.5
#
s2 line linestyle 5
s2 line color 15
s2 line linewidth 2.5
#
s3 line linestyle 4
s3 line color 4
s3 line linewidth 2.5
#
s4 line linestyle 1
s4 line color 2
s4 line linewidth 2.5
#
