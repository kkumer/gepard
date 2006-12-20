# grace batch file for fit
page size 500 350
page background fill off
arrange(2, 2, .1, .2, .1)
##############################################
with g0
##############################################
world 0, 0.08, 1.0, 45
view 0.15, 0.65, 0.6, 0.95
legend loctype view
legend 0.35, 0.92
legend length 0
legend char size 0.64
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
s0 symbol color 4
s0 symbol pattern 1
s0 symbol fill color 4
s0 symbol fill pattern 1
s0 symbol linewidth 1.0
s0 symbol linestyle 1
s0 symbol char 65
s0 symbol char font 0
s0 symbol skip 0
s0 legend  "H1, placeholderpla Q2=4"
s0 line type 0
s0 errorbar on
s0 errorbar place both
s0 errorbar color 4
s0 errorbar pattern 1
s0 errorbar size 0.49
#
s1 line linestyle 4
s1 line color 4
s1 line linewidth 1.5
#
s2 type xydy
s2 symbol 2
s2 symbol size 0.49
s2 symbol color 15
s2 symbol pattern 1
s2 symbol fill color 15
s2 symbol fill pattern 1
s2 symbol linewidth 1.0
s2 symbol linestyle 1
s2 symbol char 65
s2 symbol char font 0
s2 symbol skip 0
s2 legend  "H1, placeholderpla Q2=8"
s2 line type 0
s2 errorbar on
s2 errorbar place both
s2 errorbar color 15
s2 errorbar pattern 1
s2 errorbar size 0.49
#
s3 line linestyle 5
s3 line color 15
s3 line linewidth 1.5
#
##############################################
with g1
##############################################
world 0, 0.05, 90, 25
view 0.69, 0.65, 1.14, 0.95
legend loctype view
legend 0.84, 0.92
legend length 0
legend char size 0.64
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
s0 symbol fill pattern 1
s0 symbol linewidth 1.0
s0 symbol linestyle 1
s0 symbol char 65
s0 symbol char font 0
s0 symbol skip 0
s0 legend  "H1, placeholderpla W=82"
s0 line type 0
s0 errorbar on
s0 errorbar place both
s0 errorbar color 4
s0 errorbar pattern 1
s0 errorbar size 0.49
#
s1 line linestyle 4
s1 line color 4
s1 line linewidth 1.5
#
s2 type xydy
s2 symbol 4
s2 symbol size 0.49
s2 symbol color 2
s2 symbol pattern 1
s2 symbol fill color 2
s2 symbol fill pattern 1
s2 symbol linewidth 1.0
s2 symbol linestyle 1
s2 symbol char 65
s2 symbol char font 0
s2 symbol skip 0
s2 legend  "ZEUS, placeholderpla W=89"
s2 line type 0
s2 errorbar on
s2 errorbar place both
s2 errorbar color 15
s2 errorbar pattern 1
s2 errorbar size 0.49
#
s3 line linestyle 1
s3 line color 2
s3 line linewidth 1.5
#
##############################################
with g2
##############################################
world 30, 1., 140, 14.
view 0.15, 0.26, 0.6, 0.56
legend loctype view
legend 0.16, 0.55
legend length 0
legend char size 0.64
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
yaxes scale Normal
yaxis  label "SIGMA"
yaxis  label place spec
yaxis  label place 0.05, 0.08
yaxis  ticklabel format general
#yaxis  ticklabel prec 0
yaxis  tick major 5
yaxis  tick minor ticks 4
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
s0 errorbar color 4
s0 errorbar pattern 1
s0 errorbar size 0.49
#
s1 line linestyle 4
s1 line color 4
s1 line linewidth 1.5
#
s2 type xydy
s2 symbol 2
s2 symbol size 0.49
s2 symbol color 15
s2 symbol pattern 1
s2 symbol fill color 15
s2 symbol fill pattern 1
s2 symbol linewidth 1.0
s2 symbol linestyle 1
s2 symbol char 65
s2 symbol char font 0
s2 symbol skip 0
#s2 legend  "H1, Q2=8"
s2 line type 0
s2 errorbar on
s2 errorbar place both
s2 errorbar color 15
s2 errorbar pattern 1
s2 errorbar size 0.49
#
s3 line linestyle 5
s3 line color 15
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
s4 legend  "ZEUS, placeholderplace Q2=9.6"
s4 line type 0
s4 errorbar on
s4 errorbar place both
s4 errorbar color 2
s4 errorbar pattern 1
s4 errorbar size 0.49
#
s5 line linestyle 1
s5 line color 2
s5 line linewidth 1.5
#
##############################################
with g3
##############################################
world 9e-6, -0.1, 0.5, 0.3
view 0.69, 0.26, 1.14, 0.56
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
