# grace batch file for figNNLO
page size 500 350
page background fill off
arrange(2, 2, .1, .2, .1)
with string
string on
string char size 0.780000
string 0.43, 0.68
string def "Q225"
with string
string on
string 0.43, 0.58
string def "Q225"
with string
string on
string 0.95, 0.85
string def "Q2100"
with string
string on
string 0.95, 0.58
string def "Q2100"
##############################################
with g0
##############################################
world 9e-6, -60, 0.5, 10
view 0.15, 0.65, 0.6, 0.90
legend loctype view
legend 0.20, 0.605
legend length 6
legend char size 0.64
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel off
xaxis  tick major 10
xaxis  tick minor ticks 9
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
#  --------  y axis ------------
yaxes scale Normal
yaxis  label "DPK"
yaxis  label place spec
yaxis  label place 0.02, 0.11
yaxis  ticklabel format general
yaxis  ticklabel prec 2
yaxis  tick major 20
yaxis  tick minor ticks 3
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
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
#  --------  Lines ------------
s0 line linestyle 1
s0 line color 2
s0 line linewidth 2.5
s0 legend  "NNLO (P=2)"
#
s1 line linestyle 1
s1 line color 2
s1 line linewidth 1.5
#
s2 line linestyle 4
s2 line color 4
s2 line linewidth 2.5
s2 legend  "   NLO (P=1)"
#
s3 line linestyle 4
s3 line color 4
s3 line linewidth 1.5
#
s4 line linestyle 1
s4 line color 1
s4 line linewidth 1.0
#
##############################################
with g1
##############################################
world 9e-6, -22, 0.5, 25
view 0.69, 0.65, 1.14, 0.90
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  ticklabel off
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
#  --------  y axis ------------
yaxes scale Normal
#yaxis  label "KL1"
#yaxis  label place spec
#yaxis  label place 0.0, 0.11
yaxis  ticklabel format general
yaxis  ticklabel prec 3
yaxis  tick major 20
yaxis  tick minor ticks 3
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
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
#  --------  Lines ------------
s0 line linestyle 1
s0 line color 2
s0 line linewidth 2.5
#
s1 line linestyle 1
s1 line color 2
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
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
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
yaxis  label "DPP"
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
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
#  --------  Lines ------------
s0 line linestyle 1
s0 line color 2
s0 line linewidth 2.5
#
s1 line linestyle 1
s1 line color 2
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
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
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
yaxis  ticklabel prec 3
yaxis  tick major 0.1
yaxis  tick minor ticks 4
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
#  --------  Lines ------------
s0 line linestyle 1
s0 line color 2
s0 line linewidth 2.5
#
s1 line linestyle 1
s1 line color 2
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
