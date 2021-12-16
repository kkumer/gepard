"""Deep Inelastic Scattering."""

from math import log
import numpy as np

from . import data, mellin, model, wilson

class DISFormFactors(model.ParameterModel, mellin.MellinBarnes):
    """Deep Inelastic Scattering (DIS) Form Factors."""

    def __init__(self, **kwargs):
        # squared DIS charge factors
        if self.nf == 3:
            self.dis_charge = 2/9
        else:  # nf = 4
            self.dis_charge = 5/18
        self.wce_dis: Dict[float, np.ndarray] = {}
        super().__init__(**kwargs)


    def F2(self, pt: data.DataPoint) -> float:
        """Return DIS F2 for kinematic point."""
        try:
            wce_ar_dis = self.wce_dis[pt.Q2]
        except KeyError:
            # calculate it, first PW is the only relevant one
            wce_ar_dis = wilson.calc_wce(self, pt.Q2, 'DIS')[0, :, :]
            # memorize it for future
            self.wce_dis[pt.Q2] = wce_ar_dis
        pdf_prerot = self.gpd_H(0, 0)  # forward limit
        pdf = np.einsum('fa,ja->jf', self.frot_pdf, pdf_prerot)
        mb_int = self._dis_mellin_barnes_integral(pt.xB, wce_ar_dis, pdf)
        return self.dis_charge * mb_int / np.pi

