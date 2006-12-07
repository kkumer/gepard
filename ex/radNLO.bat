# grace batch file for figNLO
page size 500 350
page background fill off
arrange(2, 2, .1, .2, .1)
with string
string on
string char size 0.780000
string 0.47, 0.67
string def "\xD\f{}\S2\N = 0"
with string
string on
string 0.47, 0.59
string def "\xD\f{}\S2\N = 0"
with string
string on
string 0.94, 0.67
string def "\xD\f{}\S2\N = - 1 GeV\S2\N"
with string
string on
string 0.94, 0.59
string def "\xD\f{}\S2\N = - 1 GeV\S2\N"
##############################################
with g0
##############################################
world 9e-6, -80, 0.5, 0
view 0.15, 0.65, 0.6, 0.90
legend loctype view
legend 0.23, 0.605
legend length 6
legend char size 0.64
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel off
xaxis  tick major 10
xaxis  tick minor ticks 9
#  --------  y axis ------------
yaxes scale Normal
yaxis  label "D1K"
yaxis  label place spec
yaxis  label place 0.02, 0.11
yaxis  ticklabel format general
yaxis  ticklabel prec 2
yaxis  tick major 20
yaxis  tick minor ticks 3
#  --------  Lines ------------
s0 line linestyle 4
s0 line color 4
s0 line linewidth 2.5
s0 legend  "NLO, \oCS\O"
#
s1 line linestyle 4
s1 line color 4
s1 line linewidth 1.5
#
s2 line linestyle 5
s2 line color 15
s2 line linewidth 2.5
s2 legend  "NLO, \oMS\O"
#
s3 line linestyle 5
s3 line color 15
s3 line linewidth 1.5
#
##############################################
with g1
##############################################
world 9e-6, -80, 0.5, 0
view 0.69, 0.65, 1.14, 0.90
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel off
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
#  --------  y axis ------------
yaxes scale Normal
#yaxis  label "KL1"
#yaxis  label place spec
#yaxis  label place 0.0, 0.11
yaxis  ticklabel format general
yaxis  ticklabel prec 3
yaxis  tick major 20
yaxis  tick minor ticks 3
#  --------  Lines ------------
s0 line linestyle 4
s0 line color 4
s0 line linewidth 2.5
#
s1 line linestyle 4
s1 line color 4
s1 line linewidth 1.5
#
s2 line linestyle 5
s2 line color 15
s2 line linewidth 2.5
#
s3 line linestyle 5
s3 line color 15
s3 line linewidth 1.5
#
##############################################
with g2
##############################################
world 9e-6, -0.1, 0.5, 0.3
view 0.15, 0.38, 0.6, 0.63
#  --------  x axis ------------
xaxes scale Logarithmic
# xaxis  ticklabel off
xaxis  label "\xx"
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
#  --------  y axis ------------
yaxes scale Normal
yaxis  label "D1P"
yaxis  label place spec
yaxis  label place 0.0, 0.11
yaxis  ticklabel format general
yaxis  ticklabel prec 3
yaxis  label place spec
yaxis  label place 0.0, 0.11
yaxis  ticklabel format general
yaxis  ticklabel prec 3
yaxis  tick major 0.1
yaxis  tick minor ticks 4
#yaxis  tick spec type both
#yaxis  tick spec 5
#yaxis  tick major 0, 0
#yaxis  ticklabel 0, "0"
#yaxis  tick minor 1, 0.19635
#yaxis  tick major 2, 0.392699
#yaxis  ticklabel 2, "\xp\f{}/8"
#yaxis  tick minor 3, 0.589049
#yaxis  tick major 4, 0.785398
#yaxis  ticklabel 4, "\xp\f{}/4"
#  --------  Lines ------------
s0 line linestyle 4
s0 line color 4
s0 line linewidth 2.5
#
s1 line linestyle 4
s1 line color 4
s1 line linewidth 1.5
#
s2 line linestyle 5
s2 line color 15
s2 line linewidth 2.5
#
s3 line linestyle 5
s3 line color 15
s3 line linewidth 1.5
#
##############################################
with g3
##############################################
world 9e-6, -0.1, 0.5, 0.3
view 0.69, 0.38, 1.14, 0.63
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  label "\xx"
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
#  --------  y axis ------------
yaxes scale Normal
#yaxis  label "D1P"
yaxis  ticklabel format general
yaxis  ticklabel prec 3
yaxis  tick major 0.1
yaxis  tick minor ticks 4
# yaxis  tick spec type both
# yaxis  tick spec 5
# yaxis  tick major 0, 0
# yaxis  ticklabel 0, "0"
# yaxis  tick minor 1, 0.19635
# yaxis  tick major 2, 0.392699
# yaxis  ticklabel 2, "\xp\f{}/8"
# yaxis  tick minor 3, 0.589049
# yaxis  tick major 4, 0.785398
# yaxis  ticklabel 4, "\xp\f{}/4"
#  --------  Lines ------------
s0 line linestyle 4
s0 line color 4
s0 line linewidth 2.5
#
s1 line linestyle 4
s1 line color 4
s1 line linewidth 1.5
#
s2 line linestyle 5
s2 line color 15
s2 line linewidth 2.5
#
s3 line linestyle 5
s3 line color 15
s3 line linewidth 1.5
#
