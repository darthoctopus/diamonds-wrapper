import numpy as np
from os import system
from os.path import isfile, dirname

class DiamondsSession:
    def __init__(self, ID, filename, n_head=12, CATALOGUE="EPIC"):

        assert isinstance(ID, int), "Please supply a numeric-type ID number"
        
        self.ID = ID
        self.CATALOGUE = CATALOGUE

        # Set up all required files
        for _ in ['bg', 'pb']:
            system(f"mkdir -p DIAMONDS/results/{self.CATALOGUE}{ID}/{_}/")
        for _ in ['', '/pb']:
            system(f"cp DIAMONDS/config{_.replace('/', '')}/* DIAMONDS/results/{self.CATALOGUE}{ID}{_}")
        system(f"cat {filename} | tail -n +{n_head} >| DIAMONDS/data/{self.CATALOGUE}{ID}.txt")

        system(f"echo {dirname(__file__)}/DIAMONDS/ >| localPath.txt")

        self._bg = None
        self._pb = None

        self._bg_result = None
        self._pb_result = None

    @property
    def bg(self):
        '''
        MLE estimator for background parameters
        '''
        return self._bg

    @bg.setter
    def bg(self, t):
        self._bg = np.array(t)

    def write_bg_hyperparams(self, f, overwrite=False):
        assert self.bg is not None, "Need an input background model"
        fname = f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/0_bg.txt"
        if not isfile(fname) or overwrite:
            v = np.array([self.bg/f, self.bg * f]).T
            np.savetxt(fname, v, fmt="%10.10f")

    def read_bg_hyperparams(self):
        v = np.loadtxt(f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/0_bg.txt")
        return v

    def run_bg(self, force=False):
        assert isfile(f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/0_bg.txt"), "Need to have written background hyperparameters"
        outfile = f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/bg/background_parameterSummary.txt"
        if not isfile(outfile) or force:
            assert len(self.bg) == 8, "Insufficient number of parameters for TwoHarvey fit"
            system(f"DIAMONDS/bin/background {self.CATALOGUE} {self.ID} bg TwoHarvey 0 0 0 >| DIAMONDS/results/{self.CATALOGUE}{self.ID}/bg.log 2>&1")
        self._bg_result = np.loadtxt(outfile)
        self.bg = self._bg_result[:, 0]

        np.savetxt(f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/backgroundParameters.txt", self.bg[:-3], fmt="%10.10f")
        np.savetxt(f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/gaussianEnvelopeParameters.txt", self.bg[-3:], fmt="%10.10f")

    @property
    def pb(self):
        '''
        Manually identified mode frequencies (list)
        '''
        return self._pb

    @pb.setter
    def pb(self, t):
        self._pb = np.array(t)

    def write_pb_hyperparams(self, e, amplitude, width, e_amp=None, e_width=None, run=0):
        assert self.pb is not None, "Need an input linelist"
        f = np.array([self.pb - e, self.pb + e]).T
        o = np.ones_like(self.pb)
        
        if e_amp is None:
            a = np.array([o * np.sqrt(amplitude)/2, o * amplitude]).T
        else:
            a = np.array([np.maximum(0, o * (amplitude - e_amp)), o * (amplitude + e_amp)]).T

        if e_width is None:
            w = np.array([1e-5 * o * width, o * width]).T
        else:
            w = np.array([np.maximum(0, o * (width - e_width)), o * (width + e_width)]).T

        v = np.concatenate([np.array(_) for _ in zip(f, a ,w)])
        np.savetxt(f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/pb/prior_hyperParameters_{run}.txt", v, fmt="%10.10f")
        a = [np.min(self.pb- 2 * e), np.max(self.pb + 2*e)]
        np.savetxt(f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/pb/frequencyRange_{run}.txt", a, fmt="%10.10f")

    def run_pb(self, force=False, run=0):
        system(f"mkdir -p DIAMONDS/results/{self.CATALOGUE}{self.ID}/pb/{run}")
        assert isfile(f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/pb/frequencyRange_{run}.txt"), "Need to have written peakbagging hyperparameters"
        outfile = f"DIAMONDS/results/{self.CATALOGUE}{self.ID}/pb/{run}/peakbagging_parameterSummary.txt"
        if not isfile(outfile) or force:
            system(f"DIAMONDS/bin/peakbagging {self.CATALOGUE} {self.ID} pb {run} TwoHarvey prior_hyperParameters -1 0 0 0 >| DIAMONDS/results/{self.CATALOGUE}{self.ID}/pb_{run}.log 2>&1")
        self._pb_result = np.loadtxt(outfile)
        self.pb = self._pb_result[::3, 1]

