import numpy as np
import logging
from fits_reduce.util.overscancorr import OverscanCorr

import time

class T80SPreProc(OverscanCorr):

    def __init__(self,*args,**kwargs):
        OverscanCorr.__init__(self,*args,**kwargs)
        self.log = logging.getLogger("t80s_preproc")

    def gain(self):
        gain_img = np.zeros_like(self.ccd.data,dtype=np.float)

        for subarr in self._ccdsections:
            # print scan_level
            gain_img[subarr.section] += subarr.gain

        newdata = np.zeros_like(self.ccd.data,dtype=np.float) + self.ccd.data
        newdata /= gain_img
        self.ccd.data = newdata

    def norm(self):
        newdata = np.zeros_like(self.ccd.data,dtype=np.float) + self.ccd.data
        newdata /= np.median(self.ccd.data)
        self.ccd.data = newdata

    def norm_section(self):
        norm_img = np.zeros_like(self.ccd.data,dtype=np.float)

        # import pyds9 as ds9

        # d = ds9.ds9()

        level_start = np.zeros_like(self._ccdsections)
        level_end = np.zeros_like(self._ccdsections)
        n_gain = np.zeros_like(self._ccdsections) + 1.0 / np.median(self.ccd.data)

        for i in range(0,len(self._ccdsections),2):
            subarr = self._ccdsections[i]
        # for i,subarr in enumerate(self._ccdsections):
            # print scan_level
            # norm_img[subarr.section] = np.median(self.ccd.data[subarr.section])
            # d.set_np2arr(self.ccd.data[subarr.section]) #[:,-10:])
            # return
            level_start[i] = np.median(self.ccd.data[subarr.section][:,:2])
            level_end[i] = np.median(self.ccd.data[subarr.section][:,-2:])

        for i in range(1,len(self._ccdsections),2):
            subarr = self._ccdsections[i]
        # for i,subarr in enumerate(self._ccdsections):
            # print scan_level
            # norm_img[subarr.section] = np.median(self.ccd.data[subarr.section])
            # d.set_np2arr(self.ccd.data[subarr.section]) #[:,-10:])
            # return
            level_start[i] = np.median(self.ccd.data[subarr.section][:,:10])
            level_end[i] = np.median(self.ccd.data[subarr.section][:,-10:])

        for i in range(2,len(level_start),2):
            n_gain[i] = n_gain[i-2] * level_end[i-2] / level_start[i]

        for i in range(3,len(level_start),2):
            n_gain[i] = n_gain[i-2] * level_end[i-2] / level_start[i]

        print n_gain

        for i,subarr in enumerate(self._ccdsections):
            # print scan_level
            norm_img[subarr.section] = n_gain[i]

        newdata = np.zeros_like(self.ccd.data,dtype=np.float) + self.ccd.data
        newdata *= norm_img
        self.ccd.data = newdata

    def get_avg(self, filename):
        level = np.zeros((self._parallelports,self._serialports))
        for subarr in self._ccdsections:
            level[subarr.parallel_index][subarr.serial_index] = np.median(self.ccd.data[subarr.section])
        np.save(level,filename)
        return

    def linearize(self,coeffs, saturation = 60e3):
        lin_corr_img = np.zeros_like(self.ccd.data,dtype=np.float64)
        heaviside = lambda x,x0,t : 1 / (1+np.exp(-(x-x0)/t))
        def correction(x,y):
                return (x[0]+x[1]*y)*(1.-heaviside(y,x[2],x[3])) + (x[4]+x[5]*y+x[6]*y**2.+x[7]*y**3+x[8]*y**4)*heaviside(y,x[2],x[3])

        for i in range(0,len(self._ccdsections)):
            subarr = self._ccdsections[i]
            self.log.debug('Linearizing region %i x %i...' % (subarr.parallel_index,subarr.serial_index))
            if coeffs[subarr.parallel_index][subarr.serial_index]['type'] == 0:
                lin_corr_img[subarr.section] += 1.
            elif coeffs[subarr.parallel_index][subarr.serial_index]['type'] == 1:
                c = [coeffs[subarr.parallel_index][subarr.serial_index]['c1_1'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['c1_2'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['c1_3']]
                ppol = np.poly1d(c)
                lin_corr_img[subarr.section] = ppol(self.ccd.data[subarr.section])
            else:
                c = [coeffs[subarr.parallel_index][subarr.serial_index]['c1_1'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['c1_2'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['break'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['step'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['c2_1'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['c2_2'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['c2_3'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['c2_4'],
                     coeffs[subarr.parallel_index][subarr.serial_index]['c2_5']]

                lin_corr_img[subarr.section] = correction(c,self.ccd.data[subarr.section])

        mask = np.bitwise_and(lin_corr_img > 0.,
                              self.ccd.data < saturation)

        lin_corr_img[np.bitwise_not(mask)] = 1.

        newdata = self.ccd.data / lin_corr_img

        self.ccd.data = newdata

