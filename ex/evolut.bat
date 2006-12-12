# grace batch file for evolut
page size 500 350
page background fill off
arrange(2, 2, .1, .2, .1)
with string
string on
string char size 0.780000
string 0.37, 0.85
string def "Q25"
with string
string on
string 0.27, 0.41
string def "Q25"
with string
string on
string 0.78, 0.85
string def "Q2205"
with string
string on
string 0.78, 0.52
string def "Q2205"
##############################################
with g0
##############################################
world 9e-6, -10, 0.5, 60
view 0.15, 0.65, 0.6, 0.90
legend loctype view
legend 0.94, 0.87
legend length 5
legend char size 0.64
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel off
xaxis  tick major 10
xaxis  tick minor ticks 9
#  -------- alt x axis ------------
altxaxis  on
altxaxis  type zero true
altxaxis  bar on
altxaxis  bar color 1
altxaxis  bar linestyle 1
altxaxis  bar linewidth 0.5
altxaxis  tick off
altxaxis  tick major 10
altxaxis  ticklabel off
#  --------  y axis ------------
yaxes scale Normal
yaxis  label "DDK"
yaxis  label place spec
yaxis  label place 0.02, 0.11
yaxis  ticklabel format decimal
yaxis  ticklabel prec 0
yaxis  tick major 10
yaxis  tick minor ticks 4
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 2.5
s0 legend  " LO"
#
s1 line linestyle 2
s1 line color 1
s1 line linewidth 1.5
#
s2 line linestyle 4
s2 line color 4
s2 line linewidth 2.5
s2 legend  "NLO, \oCS\O"
#
s3 line linestyle 4
s3 line color 4
s3 line linewidth 1.5
#
s4 line linestyle 5
s4 line color 15
s4 line linewidth 2.5
s4 legend  "NLO, \oMS\O"
#
s5 line linestyle 5
s5 line color 15
s5 line linewidth 1.5
#
##############################################
with g1
##############################################
world 9e-6, 0, 0.5, 420
view 0.69, 0.65, 1.14, 0.90
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel off
xaxis  tick major 10
xaxis  tick minor ticks 9
#  -------- alt x axis ------------
altxaxis  on
altxaxis  type zero true
altxaxis  bar on
altxaxis  bar color 1
altxaxis  bar linestyle 1
altxaxis  bar linewidth 0.5
altxaxis  tick off
altxaxis  tick major 10
altxaxis  ticklabel off
#  --------  y axis ------------
yaxes scale Normal
#yaxis  label "KL1"
#yaxis  label place spec
#yaxis  label place 0.0, 0.11
yaxis  ticklabel format decimal
yaxis  ticklabel prec 0
yaxis  tick major 100
yaxis  tick minor ticks 4
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 2.5
#
s1 line linestyle 2
s1 line color 1
s1 line linewidth 1.5
#
s2 line linestyle 4
s2 line color 4
s2 line linewidth 2.5
#
s3 line linestyle 4
s3 line color 4
s3 line linewidth 1.5
#
s4 line linestyle 5
s4 line color 15
s4 line linewidth 2.5
#
s5 line linestyle 5
s5 line color 15
s5 line linewidth 1.5
#
##############################################
with g2
##############################################
world 9e-6, -0.1, 0.5, 0.04
view 0.15, 0.38, 0.6, 0.63
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  label "\xx"
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
#  -------- alt x axis ------------
altxaxis  on
altxaxis  type zero true
altxaxis  bar on
altxaxis  bar color 1
altxaxis  bar linestyle 1
altxaxis  bar linewidth 0.5
altxaxis  tick off
altxaxis  tick major 10
altxaxis  ticklabel off
#  --------  y axis ------------
yaxes scale Normal
yaxis  label "DDP"
yaxis  label place spec
yaxis  label place 0.0, 0.11
yaxis  ticklabel format general
yaxis  ticklabel prec 3
yaxis  label place spec
yaxis  label place 0.0, 0.11
yaxis  ticklabel format general
yaxis  ticklabel prec 2
yaxis  tick major 0.05
yaxis  tick minor ticks 4
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 2.5
#
s1 line linestyle 2
s1 line color 1
s1 line linewidth 1.5
#
s2 line linestyle 4
s2 line color 4
s2 line linewidth 2.5
#
s3 line linestyle 4
s3 line color 4
s3 line linewidth 1.5
#
s4 line linestyle 5
s4 line color 15
s4 line linewidth 2.5
#
s5 line linestyle 5
s5 line color 15
s5 line linewidth 1.5
#
##############################################
with g3
##############################################
world 9e-6, -0.28, 0.5, 0
view 0.69, 0.38, 1.14, 0.63
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  label "\xx"
xaxis  tick major 10
xaxis  tick minor ticks 9
#  -------- alt x axis ------------
altxaxis  on
altxaxis  type zero true
altxaxis  bar on
altxaxis  bar color 1
altxaxis  bar linestyle 1
altxaxis  bar linewidth 0.5
altxaxis  tick off
altxaxis  tick major 10
altxaxis  ticklabel off
#  --------  y axis ------------
yaxes scale Normal
#yaxis  label "D1P"
yaxis  ticklabel format general
yaxis  ticklabel prec 2
yaxis  tick major 0.1
yaxis  tick minor ticks 4
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 2.5
#
s1 line linestyle 2
s1 line color 1
s1 line linewidth 1.5
#
s2 line linestyle 4
s2 line color 4
s2 line linewidth 2.5
#
s3 line linestyle 4
s3 line color 4
s3 line linewidth 1.5
#
s4 line linestyle 5
s4 line color 15
s4 line linewidth 2.5
#
s5 line linestyle 5
s5 line color 15
s5 line linewidth 1.5
#
