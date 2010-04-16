# Output layer transformation function

from numpy import log

from numpy import power as Power
from numpy import sqrt as Sqrt

def BSA(ImcffH, xB, t, Q2):
    return (  (11.28*ImcffH*(1.26 - t)*(2 - (0.019307771923263693*Q2)/xB)*xB*
        (1 + (3.521417481516288*Power(xB,2))/Q2)*
        Power(1. + (3.521417481516288*Power(xB,2))/Q2,2)*
        Sqrt(-((t*(1. - 0.0003281873555975635*Q2 - 
                (0.019307771923263693*Q2)/xB)*(1. - xB)*
              (1. + (Q2*((3.521417481516288*Power(xB,2))/Q2 + 
                     2.*(1. - xB)*
                      (1. - Sqrt(1. + (3.521417481516288*Power(xB,2))/Q2))))/
                 (t*(4.*(1. - xB)*xB + (3.521417481516288*Power(xB,2))/Q2)))*
              (Sqrt(1. + (3.521417481516288*Power(xB,2))/Q2) + 
                (0.25*(4.*(1. - xB)*xB + (3.521417481516288*Power(xB,2))/Q2)*
                   (t + (Q2*((3.521417481516288*Power(xB,2))/Q2 + 
                          2.*(1. - xB)*
                           (1. - 
                             Sqrt(1. + (3.521417481516288*Power(xB,2))/Q2)))\
    )/(4.*(1. - xB)*xB + (3.521417481516288*Power(xB,2))/Q2)))/(Q2*(1. - xB)))\
    )/Q2))*(1 + ((1 - xB + (-1 + 
                  Sqrt(1 + (3.521417481516288*Power(xB,2))/Q2))/2.)*
             (t + (Q2*((3.521417481516288*Power(xB,2))/Q2 + 
                    2.*(1. - xB)*(1. - 
                       Sqrt(1. + (3.521417481516288*Power(xB,2))/Q2))))/
                (4.*(1. - xB)*xB + (3.521417481516288*Power(xB,2))/Q2)))/
           (Q2*(1 + (3.521417481516288*Power(xB,2))/Q2))))/
      (Power(0.71 - t,2)*(3.53 - t)*
        (8.*(1. - 0.0003281873555975635*Q2 - (0.019307771923263693*Q2)/xB)*
           (1. + (3.521417481516288*Power(xB,2))/Q2)*
           (-(Power(3.2/(Power(0.71 - t,2)*(3.53 - t)) + 
                  (1.41*(1.26 - t))/(Power(0.71 - t,2)*(3.53 - t)),2)*
                Power(1. - t/Q2,2)*Power(xB,2)) + 
             (7.042834963032576*(1. - 0.220088592594768*t)*
                ((1.9880999999999998*Power(1.26 - t,2))/
                   (Power(0.71 - t,4)*Power(3.53 - t,2)) - 
                  (2.9079199083179303*t)/(Power(0.71 - t,4)*Power(3.53 - t,2))\
    )*Power(xB,2))/Q2) + Power(2. - (0.019307771923263693*Q2)/xB,2)*
           (((1.9880999999999998*Power(1.26 - t,2))/
                 (Power(0.71 - t,4)*Power(3.53 - t,2)) - 
                (2.9079199083179303*t)/(Power(0.71 - t,4)*Power(3.53 - t,2)))*
              (2. + (3.521417481516288*Power(xB,2))/Q2)*
              ((3.521417481516288*Power(1. + t/Q2,2)*Power(xB,2))/t + 
                4.*(1. - xB)*(1. + (t*xB)/Q2)) + 
             4.*Power(3.2/(Power(0.71 - t,2)*(3.53 - t)) + 
                (1.41*(1.26 - t))/(Power(0.71 - t,2)*(3.53 - t)),2)*Power(xB,2)*
              (xB - (Power(t,2)*(1. - 2.*xB)*xB)/Power(Q2,2) + 
                Power(1. - t/Q2,2)*
                 (1. - xB + (1.760708740758144*Power(xB,2))/Q2))) - 
          (8.*t*(1. - 0.0003281873555975635*Q2 - (0.019307771923263693*Q2)/xB)*
             (1. - xB)*(2.*Power(3.2/(Power(0.71 - t,2)*(3.53 - t)) + 
                  (1.41*(1.26 - t))/(Power(0.71 - t,2)*(3.53 - t)),2)*
                Power(xB,2) + (Q2*
                  ((1.9880999999999998*Power(1.26 - t,2))/
                     (Power(0.71 - t,4)*Power(3.53 - t,2)) - 
                    (2.9079199083179303*t)/
                     (Power(0.71 - t,4)*Power(3.53 - t,2)))*
                  (2. + (10.564252444548865*Power(xB,2))/Q2))/t)*
             (1. + (Q2*((3.521417481516288*Power(xB,2))/Q2 + 
                    2.*(1. - xB)*(1. - 
                       Sqrt(1. + (3.521417481516288*Power(xB,2))/Q2))))/
                (t*(4.*(1. - xB)*xB + (3.521417481516288*Power(xB,2))/Q2)))*
             (Sqrt(1. + (3.521417481516288*Power(xB,2))/Q2) + 
               (0.25*(4.*(1. - xB)*xB + (3.521417481516288*Power(xB,2))/Q2)*
                  (t + (Q2*((3.521417481516288*Power(xB,2))/Q2 + 
                         2.*(1. - xB)*
                          (1. - Sqrt(1. + (3.521417481516288*Power(xB,2))/Q2))\
    ))/(4.*(1. - xB)*xB + (3.521417481516288*Power(xB,2))/Q2)))/(Q2*(1. - xB))))/
           Q2))
      )

def trans(x, theory_and_pt_and_deriv):
    theory, pt, deriv = theory_and_pt_and_deriv
    xB, t, Q2 = pt.xB, pt.t, pt.Q2
    # multiplying with (1-xB) implements constraint CFF(xB=1) = 0
    res = theory.predict(pt, parameters={'outputvalue':x})
    #res = - BSA((1-xB)*x, xB, t, Q2) 
    return res

# FIXME: this map2pt dictionary is attribute of the trans module but it would
# be natural that it is attribute of Fitter instance, related to its artificial
# data. But than it is not clear how to make it visible to brain.*
map2pt = {}
