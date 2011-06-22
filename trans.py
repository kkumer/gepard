# Output layer transformation function

def trans(o, theory_and_pt, oC=0):
    """Predict pt within theory with CFFs given by NN's output layer o.
    
    oC -- value of subtraction constant
    
    """
    theory, pt = theory_and_pt
    if theory.model.endpointpower:
        # multiplying with (1-xB)^2 to implement constraint CFF(xB=1) = 0
        o = o * (1-pt.xB)**theory.model.endpointpower
    res = theory.predict(pt, parameters={'outputvalue':o, 'outputvalueC':oC})
    return res


# FIXME: this map2pt dictionary is attribute of the trans module but it would
# be natural that it is attribute of Fitter instance, related to its artificial
# data. But than it is not clear how to make it visible to brain.*
map2pt = {}

outmem = {}
