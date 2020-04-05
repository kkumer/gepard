# Output layer transformation function

def trans(o, theory_and_pt, oC=0, oCu=0, oCd=0):
    """Predict pt within theory with CFFs given by NN's output layer o.
    
    oC -- value of subtraction constant in non-flavored model
    oCu, oCd -- values of subtraction constant in flavored model
    
    """
    theory, pt = theory_and_pt
    if theory.model.endpointpower:
        # multiplying imag parts of CFF with (1-xB)^2 to make them zero at xB=1
        fac = (1-pt.xB)**theory.model.endpointpower
        ep = [int(cff[:2]=='Re')+int(cff[:2]=='Im')*fac for cff in theory.m.output_layer]
        o = o * ep
    if theory.model.zeropointpower:
        # multiplying Eb with 1/xi if Eb and not Et is parametrized by NN
        zp = [int(cff[2:]!='Et')+int(cff[2:]=='Et')*(2.-pt.xB)/pt.xB 
                for cff in theory.m.output_layer]
        o = o * zp
    res = theory.predict(pt, parameters={'outputvalue':o, 'outputvalueC':oC, 
        'outputvalueCu':oCu, 'outputvalueCd':oCd})
    return res


# FIXME: this map2pt dictionary is attribute of the trans module but it would
# be natural that it is attribute of Fitter instance, related to its artificial
# data. But than it is not clear how to make it visible to brain.*
map2pt = {}

outmem = {}
