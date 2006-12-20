# grace batch file for figQ
page size 500 350
page background fill off
arrange(1, 2, .1, .2, .1)
with string
string on
string char size 0.90
string 0.18, 0.81
string def "xi5"
with string
string on
string char size 0.90
string 0.35, 0.725
string def "xi3"
with string
string on
string char size 0.90
string 0.35, 0.623
string def "xi2"
with string
string on
string char size 0.90
string 0.35, 0.514
string def "xi1"
with string
string on
string char size 0.90
string 0.95, 0.511
string def "xi1a"
with string
string on
string char size 0.90
string 0.95, 0.593
string def "xi2a"
with string
string on
string char size 0.90
string 0.95, 0.672
string def "xi3a"
with string
string on
string char size 0.90
string 0.95, 0.786
string def "xi5a"
##############################################
with g0
##############################################
world 0, 3, 100, 110
view 0.15, 0.45, 0.6, 0.85
# legend loctype view
# legend 0.23, 0.605
# legend length 7
# legend char size 0.64
#  --------  x axis ------------
xaxes scale Normal
xaxis  label "calQ2"
xaxis  ticklabel format general
xaxis  tick major 20
xaxis  tick minor ticks 3
#  --------  y axis ------------
yaxes scale Logarithmic 
yaxis  label "XIH"
yaxis  label place spec
yaxis  label place 0.0, 0.1
yaxis  tick spec type both
yaxis  tick spec 9
yaxis  tick major 0, 5
yaxis  ticklabel 0, "5"
yaxis  tick minor 1, 7.5
yaxis  tick major 2, 10
yaxis  ticklabel 2, "10"
yaxis  tick minor 3, 15
yaxis  tick major 4, 20
yaxis  ticklabel 4, "20"
yaxis  tick minor 5, 35
yaxis  tick major 6, 50
yaxis  ticklabel 6, "50"
yaxis  tick minor 7, 75
yaxis  tick major 8, 100
yaxis  ticklabel 8, "100"
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 1.5
# s0 legend  "  LO"
#
s1 line linestyle 4
s1 line color 4
s1 line linewidth 1.5
# s1 legend  " NLO"
#
s2 line linestyle 1
s2 line color 2
s2 line linewidth 1.5
# s1 legend  "NNLO"
#
s3 line linestyle 2
s3 line color 1
s3 line linewidth 1.5
#
s4 line linestyle 4
s4 line color 4
s4 line linewidth 1.5
#
s5 line linestyle 1
s5 line color 2
s5 line linewidth 1.5
#
s6 line linestyle 2
s6 line color 1
s6 line linewidth 1.5
#
s7 line linestyle 4
s7 line color 4
s7 line linewidth 1.5
#
s8 line linestyle 1
s8 line color 2
s8 line linewidth 1.5
#
s9 line linestyle 2
s9 line color 1
s9 line linewidth 1.5
#
s10 line linestyle 4
s10 line color 4
s10 line linewidth 1.5
#
s11 line linestyle 1
s11 line color 2
s11 line linewidth 1.5
#
##############################################
with g1
##############################################
world 0, 0.9, 100, 1.8
view 0.73, 0.45, 1.18, 0.85
#  --------  x axis ------------
xaxes scale Normal
xaxis  label "calQ2"
xaxis  ticklabel format general
xaxis  ticklabel prec 0
xaxis  tick major 20
xaxis  tick minor ticks 9
#  --------  y axis ------------
yaxes scale Normal
yaxis  label "ARGH"
yaxis  label place spec
yaxis  label place 0.0, 0.1
yaxis  ticklabel format general
yaxis  ticklabel prec 3
yaxis  tick major 0.2
yaxis  tick minor ticks 3
#  --------  Lines ------------
s0 line linestyle 2
s0 line color 1
s0 line linewidth 1.5
# s0 legend  "  LO"
#
s1 line linestyle 4
s1 line color 4
s1 line linewidth 1.5
# s1 legend  " NLO"
#
s2 line linestyle 1
s2 line color 2
s2 line linewidth 1.5
# s1 legend  "NNLO"
#
s3 line linestyle 2
s3 line color 1
s3 line linewidth 1.5
#
s4 line linestyle 4
s4 line color 4
s4 line linewidth 1.5
#
s5 line linestyle 1
s5 line color 2
s5 line linewidth 1.5
#
s6 line linestyle 2
s6 line color 1
s6 line linewidth 1.5
#
s7 line linestyle 4
s7 line color 4
s7 line linewidth 1.5
#
s8 line linestyle 1
s8 line color 2
s8 line linewidth 1.5
#
s9 line linestyle 2
s9 line color 1
s9 line linewidth 1.5
#
s10 line linestyle 4
s10 line color 4
s10 line linewidth 1.5
#
s11 line linestyle 1
s11 line color 2
s11 line linewidth 1.5
#
