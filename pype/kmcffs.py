import shelve, Model, Data

db = shelve.open('theories.db',  flag='r')
th = db['KM15']

cffs = ['ImH', 'ReH', 'ImE', 'ReE', 'ImHt', 'ReHt', 'ImEt', 'ReEt']


pt = Data.DummyPoint()

def kmcffs(xB, t, Q2):
    pt.Q2 = Q2
    pt.t = t
    pt.xB = xB
    pt.xi = xB/(2.-xB)

    row = [pt.xB, pt.t, pt.Q2]
    for cff in cffs:
        th.m.g.newcall = 1
        row.append(getattr(th.m, cff)(pt))

    return row


#print(getcffs(0.36, -0.2))
