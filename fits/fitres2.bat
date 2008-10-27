# grace batch file for fitres2
page size 500 350
page background fill off
arrange(2, 2, .1, .2, .1)
##############################################
with g0
##############################################
world 0, 0.03, 1.0, 20
view 0.15, 0.65, 0.6, 0.95
legend loctype view
legend 0.33, 0.93
legend length 0
legend char size 0.64
legend box linestyle 0
legend box fill pattern 0
#  --------  x axis ------------
xaxes scale Normal
xaxis  label "MT"
xaxis  ticklabel prec 1
xaxis  tick major 0.2
xaxis  tick minor ticks 3
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
#  --------  y axis ------------
yaxes scale Logarithmic
yaxis  label "PARSIGMA"
yaxis  label place spec
yaxis  label place 0.07, 0.08
yaxis  ticklabel format general
#yaxis  ticklabel prec 2
yaxis  tick major 10
yaxis  tick minor ticks 9
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
#  --------  Lines ------------
s0 type xydy
s0 symbol 1
s0 symbol size 0.49
s0 symbol color 1
s0 symbol pattern 1
s0 symbol fill color 1
s0 symbol fill pattern 1
s0 symbol linewidth 1.0
s0 symbol linestyle 1
s0 symbol char 65
s0 symbol char font 0
s0 symbol skip 0
s0 legend  "H1, placeholderpla Q2=8"
s0 line type 0
s0 errorbar on
s0 errorbar place both
s0 errorbar color 1
s0 errorbar pattern 1
s0 errorbar size 0.49
#
s1 line linestyle 1
s1 line color 1
s1 line linewidth 1.5
#
s2 type xydy
s2 symbol 4
s2 symbol size 0.49
s2 symbol color 4
s2 symbol pattern 1
s2 symbol fill color 4
s2 symbol fill pattern 1
s2 symbol linewidth 1.0
s2 symbol linestyle 1
s2 symbol char 65
s2 symbol char font 0
s2 symbol skip 0
s2 legend  "H1, placeholderpla Q2=15.5"
s2 line type 0
s2 errorbar on
s2 errorbar place both
s2 errorbar color 4
s2 errorbar pattern 1
s2 errorbar size 0.49
#
s3 line linestyle 1
s3 line color 4
s3 line linewidth 1.5
#
s4 type xydy
s4 symbol 2
s4 symbol size 0.49
s4 symbol color 12
s4 symbol pattern 1
s4 symbol fill color 12
s4 symbol fill pattern 1
s4 symbol linewidth 1.0
s4 symbol linestyle 1
s4 symbol char 65
s4 symbol char font 0
s4 symbol skip 0
s4 legend  "H1, placeholderpla Q2=25"
s4 line type 0
s4 errorbar on
s4 errorbar place both
s4 errorbar color 12
s4 errorbar pattern 1
s4 errorbar size 0.49
#
s5 line linestyle 1
s5 line color 12
s5 line linewidth 1.5
#
##############################################
with g1
##############################################
world 0, 0.05, 90, 25
view 0.69, 0.65, 1.14, 0.95
legend loctype view
legend 0.79, 0.93
legend length 0
legend char size 0.64
legend box linestyle 0
legend box fill pattern 0
#  --------  x axis ------------
xaxes scale Normal
xaxis  label "Q2"
xaxis  ticklabel format general
#xaxis  ticklabel prec 0
xaxis  tick major 20
xaxis  tick minor ticks 3
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
#  --------  y axis ------------
yaxes scale Logarithmic
yaxis  label "SIGMA"
yaxis  label place spec
yaxis  label place 0.05, 0.07
yaxis  ticklabel format general
yaxis  tick major 10
yaxis  tick minor ticks 9
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
#  --------  Lines ------------
s0 type xydy
s0 symbol 1
s0 symbol size 0.49
s0 symbol color 4
s0 symbol pattern 1
s0 symbol fill color 4
s0 symbol fill pattern 0
s0 symbol linewidth 1.0
s0 symbol linestyle 1
s0 symbol char 65
s0 symbol char font 0
s0 symbol skip 0
s0 legend  "H1, HI, placeholder W=82"
s0 line type 0
s0 errorbar on
s0 errorbar place both
s0 errorbar color 4
s0 errorbar pattern 1
s0 errorbar size 0.49
#
s1 line linestyle 1
s1 line color 4
s1 line linewidth 1.5
#
s2 type xydy
s2 symbol 1
s2 symbol size 0.49
s2 symbol color 4
s2 symbol pattern 1
s2 symbol fill color 4
s2 symbol fill pattern 1
s2 symbol linewidth 1.0
s2 symbol linestyle 1
s2 symbol char 65
s2 symbol char font 0
s2 symbol skip 0
s2 legend  "H1, HII, placeholder W=82"
s2 line type 0
s2 errorbar on
s2 errorbar place both
s2 errorbar color 4
s2 errorbar pattern 1
s2 errorbar size 0.49
#
s3 line linestyle 1
s3 line color 4
s3 line linewidth 1.5
#
s4 type xydy
s4 symbol 4
s4 symbol size 0.49
s4 symbol color 2
s4 symbol pattern 1
s4 symbol fill color 2
s4 symbol fill pattern 1
s4 symbol linewidth 1.0
s4 symbol linestyle 1
s4 symbol char 65
s4 symbol char font 0
s4 symbol skip 0
s4 legend  "ZEUS, placeholderpla W=89"
s4 line type 0
s4 errorbar on
s4 errorbar place both
s4 errorbar color 15
s4 errorbar pattern 1
s4 errorbar size 0.49
#
s5 line linestyle 5
s5 line color 2
s5 line linewidth 1.5
#
##############################################
with g2
##############################################
world 30, 0.1, 140, 10.
view 0.15, 0.26, 0.6, 0.56
legend loctype view
legend 0.16, 0.55
legend length 0
legend char size 0.64
legend box linestyle 0
legend box fill pattern 0
legend off
#  --------  x axis ------------
xaxes scale Normal
# xaxis  ticklabel off
xaxis  label "WW"
xaxis  ticklabel format general
#xaxis  ticklabel prec 0
xaxis  tick major 40
xaxis  tick minor ticks 3
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
#  --------  y axis ------------
yaxes scale Logarithmic
yaxis  label "SIGMA"
yaxis  label place spec
yaxis  label place 0.05, 0.08
yaxis  ticklabel format general
#yaxis  ticklabel prec 0
yaxis  tick major 10
yaxis  tick minor ticks 9
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
#  --------  Lines ------------
s0 type xydy
s0 symbol 1
s0 symbol size 0.49
s0 symbol color 1
s0 symbol pattern 1
s0 symbol fill color 1
s0 symbol fill pattern 1
s0 symbol linewidth 1.0
s0 symbol linestyle 1
s0 symbol char 65
s0 symbol char font 0
s0 symbol skip 0
s0 legend  "H1, Q2=8"
s0 line type 0
s0 errorbar on
s0 errorbar place both
s0 errorbar color 1
s0 errorbar pattern 1
s0 errorbar size 0.49
#
s1 line linestyle 1
s1 line color 1
s1 line linewidth 1.5
#
s2 type xydy
s2 symbol 4
s2 symbol size 0.49
s2 symbol color 4
s2 symbol pattern 1
s2 symbol fill color 4
s2 symbol fill pattern 1
s2 symbol linewidth 1.0
s2 symbol linestyle 1
s2 symbol char 65
s2 symbol char font 0
s2 symbol skip 0
s2 legend  "H1, Q2=15.5"
s2 line type 0
s2 errorbar on
s2 errorbar place both
s2 errorbar color 4
s2 errorbar pattern 1
s2 errorbar size 0.49
#
s3 line linestyle 1
s3 line color 4
s3 line linewidth 1.5
#
s4 type xydy
s4 symbol 2
s4 symbol size 0.49
s4 symbol color 12
s4 symbol pattern 1
s4 symbol fill color 12
s4 symbol fill pattern 1
s4 symbol linewidth 1.0
s4 symbol linestyle 1
s4 symbol char 65
s4 symbol char font 0
s4 symbol skip 0
s4 legend  "H1, Q2=25"
s4 line type 0
s4 errorbar on
s4 errorbar place both
s4 errorbar color 12
s4 errorbar pattern 1
s4 errorbar size 0.49
#
s5 line linestyle 1
s5 line color 12
s5 line linewidth 1.5
#
##############################################
with g3
##############################################
world 3, 0.5, 94, 1.6
view 0.69, 0.26, 1.14, 0.56
#  --------  x axis ------------
xaxes scale Logarithmic
xaxis  label "Q2"
xaxis  ticklabel format power
xaxis  ticklabel prec 0
xaxis  tick major 10
xaxis  tick minor ticks 9
xaxis  tick major size 0.6
xaxis  tick minor size 0.3
xaxis  ticklabel char size 0.86
xaxis  tick spec type both
xaxis  tick spec 10
xaxis  tick major 0, 5
xaxis  ticklabel 0, "5"
xaxis  tick major 1, 10
xaxis  ticklabel 1, "10"
xaxis  tick major 2, 20
xaxis  ticklabel 2, "20"
xaxis  tick major 3, 50
xaxis  ticklabel 3, "50"
xaxis  tick major 4, 90
xaxis  ticklabel 4, "90"
xaxis  tick minor 5, 80
xaxis  tick minor 6, 70
xaxis  tick minor 7, 60
xaxis  tick minor 8, 40
xaxis  tick minor 9, 30
#  --------  y axis ------------
yaxes scale Logarithmic
yaxis  label "F2"
yaxis  ticklabel format decimal
yaxis  ticklabel prec 0
yaxis  tick major 20
yaxis  tick minor ticks 4
yaxis  tick major size 0.6
yaxis  tick minor size 0.3
yaxis  ticklabel char size 0.86
yaxis  tick spec type both
yaxis  tick spec 11
yaxis  tick major 0, 0.6
yaxis  ticklabel 0, "0.6"
yaxis  tick major 1, 1.0
yaxis  ticklabel 1, "1.0"
yaxis  tick major 2, 1.4
yaxis  ticklabel 2, "1.4"
yaxis  tick minor 3, 0.7
yaxis  tick minor 4, 0.8
yaxis  tick minor 5, 0.9
yaxis  tick minor 6, 1.1
yaxis  tick minor 7, 1.2
yaxis  tick minor 8, 1.3
yaxis  tick minor 9, 1.5
yaxis  tick minor 10, 0.5
#  --------  Lines ------------
s0 type xydy
s0 symbol 3
s0 symbol size 0.49
s0 symbol color 13
s0 symbol pattern 1
s0 symbol fill color 13
s0 symbol fill pattern 1
s0 symbol linewidth 1.0
s0 symbol linestyle 1
s0 symbol char 65
s0 symbol char font 0
s0 symbol skip 0
#s0 legend  "H1, Q2=4"
s0 line type 0
s0 errorbar on
s0 errorbar place both
s0 errorbar color 13
s0 errorbar pattern 1
s0 errorbar size 0.49
#
s1 line linestyle 1
s1 line color 13
s1 line linewidth 1.5
#
s2 type xydy
s2 symbol 7
s2 symbol size 0.49
s2 symbol color 8
s2 symbol pattern 1
s2 symbol fill color 8
s2 symbol fill pattern 1
s2 symbol linewidth 1.0
s2 symbol linestyle 1
s2 symbol char 65
s2 symbol char font 0
s2 symbol skip 0
#s2 legend  "H1, Q2=4"
s2 line type 0
s2 errorbar on
s2 errorbar place both
s2 errorbar color 8
s2 errorbar pattern 1
s2 errorbar size 0.49
#
s3 line linestyle 1
s3 line color 8
s3 line linewidth 1.5
#
s4 type xydy
s4 symbol 3
s4 symbol size 0.49
s4 symbol color 13
s4 symbol pattern 1
s4 symbol fill color 13
s4 symbol fill pattern 1
s4 symbol linewidth 1.0
s4 symbol linestyle 1
s4 symbol char 65
s4 symbol char font 0
s4 symbol skip 0
#s4 legend  "H1, Q2=4"
s4 line type 0
s4 errorbar on
s4 errorbar place both
s4 errorbar color 13
s4 errorbar pattern 1
s4 errorbar size 0.49
#
s5 line linestyle 1
s5 line color 13
s5 line linewidth 1.5
#
s6 type xydy
s6 symbol 7
s6 symbol size 0.49
s6 symbol color 8
s6 symbol pattern 1
s6 symbol fill color 8
s6 symbol fill pattern 1
s6 symbol linewidth 1.0
s6 symbol linestyle 1
s6 symbol char 65
s6 symbol char font 0
s6 symbol skip 0
#s6 legend  "H1, Q2=4"
s6 line type 0
s6 errorbar on
s6 errorbar place both
s6 errorbar color 8
s6 errorbar pattern 1
s6 errorbar size 0.49
#
s7 line linestyle 1
s7 line color 8
s7 line linewidth 1.5
#
s8 type xydy
s8 symbol 2
s8 symbol size 0.49
s8 symbol color 15
s8 symbol pattern 1
s8 symbol fill color 15
s8 symbol fill pattern 1
s8 symbol linewidth 1.0
s8 symbol linestyle 1
s8 symbol char 65
s8 symbol char font 0
s8 symbol skip 0
#s8 legend  "H1, Q2=4"
s8 line type 0
s8 errorbar on
s8 errorbar place both
s8 errorbar color 15
s8 errorbar pattern 1
s8 errorbar size 0.49
#
s9 line linestyle 1
s9 line color 15
s9 line linewidth 1.5
#
s10 type xydy
s10 symbol 1
s10 symbol size 0.49
s10 symbol color 4
s10 symbol pattern 1
s10 symbol fill color 4
s10 symbol fill pattern 1
s10 symbol linewidth 1.0
s10 symbol linestyle 1
s10 symbol char 65
s10 symbol char font 0
s10 symbol skip 0
#s10 legend  "H1, Q2=4"
s10 line type 0
s10 errorbar on
s10 errorbar place both
s10 errorbar color 4
s10 errorbar pattern 1
s10 errorbar size 0.49
#
s11 line linestyle 4
s11 line color 4
s11 line linewidth 1.5
#
##############################################
#  strings                                   #
##############################################
with string
string on
string 0.56, 0.66
string char size 1.0
string def "(a)"
with string
string on
string 1.10, 0.76
string char size 1.0
string def "(b)"
with string
string on
string 0.56, 0.28
string char size 1.0
string def "(c)"
with string
string on
string 1.09, 0.28
string char size 1.0
string def "(d)"
with string
string on
string 0.18, 0.68
string color 1
string char size 0.7
string def "W = 82 GeV"
with string
string on
string 0.27, 0.28
string color 1
string char size 0.7
string def "legend: same as on (a)"
