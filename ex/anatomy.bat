# grace batch file for figNLO
page size 500 350
page background fill off
arrange(2, 2, .1, .2, .1)
with string
string on
string char size 0.780000
string 0.31, 0.67
string def "soft"
with string
string on
string 0.31, 0.59
string def "soft"
with string
string on
string 0.94, 0.67
string def "hard"
with string
string on
string 0.94, 0.59
string def "hard"
##############################################
with g0
##############################################
world 9e-6, 8e-3, 0.5, 3.0
view 0.15, 0.65, 0.6, 0.90
legend loctype view
legend 0.37, 0.59
legend length 6
legend char size 0.64
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel off
xaxis  tick major 10
xaxis  tick minor ticks 9
#  --------  y axis ------------
yaxes scale Logarithmic
yaxis  label "DMOD"
yaxis  label place spec
yaxis  label place 0.02, 0.11
yaxis  ticklabel format power
yaxis  ticklabel prec 0
yaxis  tick major 10
yaxis  tick minor ticks 9
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 1.5
s0 legend  "DELWcoe"
#
s1 line linestyle 4
s1 line color 1
s1 line linewidth 1.5
s1 legend  "DELevoD"
#
s2 line linestyle 1
s2 line color 2
s2 line linewidth 1.5
s2 legend  "DELevoND"
#
##############################################
with g1
##############################################
world 9e-6, 8e-3, 0.5, 3.0
view 0.69, 0.65, 1.14, 0.90
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel off
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
#  --------  y axis ------------
yaxes scale Logarithmic
# yaxis  label "D1K"
# yaxis  label place spec
# yaxis  label place 0.02, 0.11
yaxis  ticklabel format power
yaxis  ticklabel prec 0
yaxis  tick major 10
yaxis  tick minor ticks 9
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 2.5
#
s1 line linestyle 4
s1 line color 1
s1 line linewidth 2.5
#
s2 line linestyle 1
s2 line color 2
s2 line linewidth 2.5
#
##############################################
with g2
##############################################
world 9e-6, -0.2, 0.5, 0.5
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
yaxis  label "DARG"
yaxis  label place spec
yaxis  label place 0.0, 0.11
yaxis  ticklabel format general
yaxis  ticklabel prec 1
yaxis  tick major 0.2
yaxis  tick minor ticks 3
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 1.5
#
s1 line linestyle 4
s1 line color 1
s1 line linewidth 1.5
#
s2 line linestyle 1
s2 line color 2
s2 line linewidth 1.5
#
##############################################
with g3
##############################################
world 9e-6, -0.2, 0.5, 0.5
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
yaxis  ticklabel prec 1
yaxis  tick major 0.2
yaxis  tick minor ticks 3
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 2.5
#
s1 line linestyle 4
s1 line color 1
s1 line linewidth 2.5
#
s2 line linestyle 1
s2 line color 2
s2 line linewidth 2.5
#
