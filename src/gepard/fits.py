"""Fit results."""

import numpy as np

from . import cff, data, dvcs, eff, gpd

GLOpoints = data.dset[31][12:] + data.dset[8] + data.dset[29]

DVCSpoints = data.dset[36] + data.dset[37] + data.dset[38] + data.dset[39] + \
  data.dset[40] + data.dset[41] + data.dset[42] + data.dset[43] + \
  data.dset[44] + data.dset[45]

H1ZEUSpoints = DVCSpoints + data.dset[48]

ALTGLO5points = data.dset[5] + data.dset[8] + data.dset[32][18:]   # DM's CLAS ALU
XLUwpoints = data.select(data.dset[50], criteria=['FTn == -1'])
XUUwpoints = data.select(data.dset[51], criteria=['FTn>=0', 'FTn <= 1'])
UNP5points = ALTGLO5points + XUUwpoints + XLUwpoints

H1ZEUSindependentNEW = data.dset[45] + data.dset[39] + data.dset[63] + data.dset[46]
H1ZEUS = H1ZEUSindependentNEW + data.select(data.dset[47], criteria=['Q2 >= 4.0'])

ALUIpoints = data.select(data.dset[68], criteria=['FTn == -1'])  # HERMES
AC0points = data.select(data.dset[67], criteria=['FTn == 0'])  # HERMES
AC1points = data.select(data.dset[67], criteria=['FTn == 1'])  # HERMES
ALUIpts = ALUIpoints[:6]
ACpts = AC0points[:6] + AC1points[:6]

TSA1points = data.select(data.dset[52], criteria=['FTn == -1'])  # HERMES A_UL
H_AULpts = TSA1points[:4]
C_AULpts = data.dset[54][:3]
AULpts = H_AULpts + C_AULpts

BTSApoints = data.select(data.dset[53], criteria=['FTn==0'])   # HERMES A_LL
ALLpts = BTSApoints[:4]

AUTIpoints = data.select(data.dset[66], criteria=['FTn==1'])  # \sin\varphi\cos\phi
AUTIpts = AUTIpoints[:4]

CLAS08pts = data.select(data.dset[81], criteria=['FTn == -1'])[-3:]

GLO15new = data.dset[94]+data.dset[95]+data.dset[96]+data.dset[101] + \
           data.dset[102]+data.dset[116]+data.dset[117]

GLO15b = H1ZEUS + ALUIpts + ACpts + CLAS08pts + AULpts + ALLpts + AUTIpts + GLO15new


class KM09(eff.DipoleEFF, cff.DispersionFixedPoleCFF, dvcs.hotfixedBMK):  # noqa: D101
    pass


th_KM09a = KM09()
par_KM09a = {'tMv': 0.8, 'rS': 1.0, 'alv': 0.43, 'Nv': 1.35,
             'rv': 0.9479, 'NS': 1.5, 'alS': 1.13, 'alpS': 0.15,
             'C': 0.24, 'tNv': 0.0, 'tbv': 3.0, 'bv': 0.4481, 'mv2': 0.64, 'alpv': 0.85,
             'mC2': 0.25020004, 'mS2': 0.5, 'bS': 3.0879, 'trv': 0.0}
th_KM09a.parameters.update(par_KM09a)
th_KM09a.name = 'KM09a'
pts_KM09a = GLOpoints

th_KM09b = KM09()
par_KM09b = {'tmv2': 0.64, 'rS': 1.0, 'alv': 0.43, 'Nv': 1.35, 'rv': 1.1064, 'NS': 1.5,
             'alS': 1.13, 'alpS': 0.15, 'C': 6.042, 'tNv': 0.6, 'tbv': 1.5, 'bv': 2.398,
             'mv2': 0.64, 'alpv': 0.85, 'mC2': 2.32196644, 'mS2': 0.5,
             'bS': 4.5845, 'trv': 4.8966}
th_KM09b.parameters.update(par_KM09b)
th_KM09b.name = 'KM09b'
pts_KM09b = GLOpoints + data.dset[30]


class KM10(eff.DipoleEFF, gpd.PWNormGPD, cff.HybridFreePoleCFF, dvcs.BM10):  # noqa: D101, E501
    pass


th_KM10 = KM10()
par_KM10 = {'tmv2': 0.7812091281581667, 'rS': 1.0, 'rpi': 3.5355742996824659,
            'alv': 0.43, 'mpi2': 0.528149679167846, 'Nv': 1.35,
            'rv': 0.61981779615528687, 'Nsea': 0.0, 'alS': 1.13, 'alpS': 0.15,
            'C': 8.7771705279055432, 'tNv': 0.6, 'tbv': 2.0496515976221952,
            'bv': 0.40375381018570683, 'mv2': 15.999972534200204, 'alpv': 0.85,
            'mC2': 0.949933596064363, 'bS': 2.0, 'trv': 7.7592155354071064,
            'delms2': 0.0, 'delmg2': 0.0,
            'ns':  0.15203911208796006, 'al0s': 1.1575060246398083, 'alps': 0.15,
            'ms2': 0.51318802297356658, 'secs': 0.27800486286895021, 'ng': 0.5,
            'this': -0.12999801345190121,
            'al0g': 1.247316701070471, 'alpg': 0.15, 'mg2': 0.7,
            'secg': -2.787206139384161, 'thig': 0.89633675296304127,
            'kaps': 0.0, 'kapg': 0.0}
th_KM10.parameters.update(par_KM10)
th_KM10.name = 'KM10'
pts_KM10 = H1ZEUSpoints + UNP5points


class KM10b(eff.KellyEFF, gpd.PWNormGPD, cff.HybridFreePoleCFF, dvcs.BM10):  # noqa: D101, E501
    pass


th_KM10b = KM10b()
par_KM10b = {'tmv2': 0.64, 'rS': 1.0, 'rpi': 4.0201, 'alv': 0.43, 'Nsea': 0.0,
             'Nv': 1.35, 'rv': 0.8081, 'mpi2': 2.36206161, 'alS': 1.13, 'alpS': 0.15,
             'C': 5.4259, 'tNv': 0.6, 'bS': 2.0, 'bv': 0.7706, 'mv2': 0.64,
             'tbv': 1.0, 'alpv': 0.85, 'mC2': 1.77023025, 'trv': 3.2931,
             'ns': 0.15203911208796006,
             'al0s': 1.1575060246398083, 'alps': 0.15, 'al0g': 1.247316701070471,
             'secs': -0.4600511871918772, 'this': 0.09351798951979662,
             'alpg': 0.15, 'mg2': 0.7, 'secg': -2.5151319493485427,
             'ms2': 0.49754317018981614,
             'thig': 0.8915757559175185, 'kaps': 0.0, 'kapg': 0.0}
th_KM10b.parameters.update(par_KM10b)
th_KM10b.name = 'KM10b'
pts_KM10b = DVCSpoints+GLOpoints+data.dset[30]


class AFKM12(eff.KellyEFF, gpd.PWNormGPD, cff.MellinBarnesCFF, dvcs.BM10):  # noqa: D101

    def E(self, eta: float, t: float) -> np.ndarray:
        """Return (npts, 4) array E_j^a for all j-points and 4 flavors."""
        # Implement BS+BG=0 sum rule that fixes 'kapg'
        self.parameters['kapg'] = - self.parameters['kaps'] * self.parameters['ns'] / (
                0.6 - self.parameters['ns'])
        kappa = np.array([self.parameters['kaps'], self.parameters['kapg'], 0, 0])
        return kappa * gpd.singlet_ng_constrained_E(
                self.jpoints, t, self.parameters, self.residualt).transpose()


th_AFKM12 = AFKM12(residualt='exp')
# In AFKM12 paper we use b instead of m2, where m2 = 1/(2b)
par_AFKM12 = {"ns": 0.152, "al0s": 1.1575, "alps": 0.1, "ms2": 0.17857142857142858,
              "delms2": 0,
              "pows": 2., "secs": 0.5127, "this": -0.21,
              "al0g": 1.2473, "alpg": 0.1, "mg2": 0.25, "delmg2": 0,
              "powg": 2., "secg": -4.8055, "thig": 1.8638,
              "kaps": 1.5,
              "Ens": 0.152,  # for fitting should be pegged to ns
              "Eal0s": 1.1575, "Ealps": 0.02, "Ems2": 0.17857142857142858, "Edelms2": 0,
              "Epows": 2., "Esecs": 0.5127, "Ethis": -0.21,
              "Eal0g": 1.2473, "Ealpg": 0.05, "Emg2": 0.25, "Edelmg2": 0,
              "Epowg": 2., "Esecg": -4.8055, "Ethig": 1.8638}
th_AFKM12.parameters.update(par_AFKM12)
th_AFKM12.name = 'AFKM12'
pts_AFKM12 = H1ZEUS


class KM15(eff.KellyEFF, gpd.PWNormGPD, cff.HybridFreePoleCFF, dvcs.BM10tw2):  # noqa: D101, E501
    pass


th_KM15 = KM15()
par_KM15 = {'tmv2': 15.94293227053628, 'rS': 1.0, 'alv': 0.43, 'tal': 0.43,
            'mpi2': 15.999998816676918, 'Nv': 1.35, 'rv': 0.918393047884448,
            'Nsea': 0.0, 'alS': 1.13, 'rpi': 2.6463144464701536,  'alpS': 0.15,
            'C': 2.7678681812890016, 'tNv': 0.6, 'bS': 2.0, 'tbv': 0.4000000003259146,
            'bv': 0.40000206775282354, 'mv2': 0.6228490959696493, 'alpv': 0.85,
            'talp': 0.85, 'mC2': 1.4497411858248308, 'trv': 0.881085721967267,
            'pows': 2.0, 'Ealpg': 0.15, 'md2': 1.0, 'mg2': 0.7, 'Ealps': 0.15,
            'powg': 2.0,
            'al0d': 0.5, 'delms2': 0.0,  'delmg2': 0.0, 'nd': 1.0, 'alpd': 1.0,
            'ns':  0.15203911208796006, 'al0s': 1.1575060246398083, 'alps': 0.15,
            'ms2': 0.4818827240886959, 'secs': 1.0707825621025808, 'ng': 0.5,
            'this': -0.36618269477432946, 'al0g': 1.247316701070471, 'alpg': 0.15,
            'mg2': 0.7, 'secg': -2.990809378821039, 'thig': 0.9052207712570559,
            'kaps': 0.0, 'kapg': 0.0}
th_KM15.parameters.update(par_KM15)
# Fitting parameters:
th_KM15._release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi',
                            'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')
th_KM15.name = 'KM15'
pts_KM15 = GLO15b
