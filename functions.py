import warnings
import numpy as np
import math

warnings.filterwarnings("ignore", category=RuntimeWarning)


def vector_coding(proximal, distal):
    thetaP = proximal
    thetaD = distal
    gamma = np.zeros((thetaP.shape[0]-1, thetaP.shape[1]))
    for indexStride in range(thetaP.shape[1]):
        for indexInstant in range(thetaP.shape[0]-1):
            thetaPDiff = thetaP[indexInstant+1, indexStride] - thetaP[indexInstant, indexStride]
            thetaDDiff = thetaD[indexInstant+1, indexStride] - thetaD[indexInstant, indexStride]
            if thetaPDiff > 0:
                gamma[indexInstant, indexStride] = math.degrees(math.atan(thetaDDiff/thetaPDiff))
            else:
                gamma[indexInstant, indexStride] = math.degrees(math.atan(thetaDDiff/thetaPDiff))+180
            if (thetaPDiff == 0) and (thetaDDiff > 0):
                gamma[indexInstant, indexStride] = 90
            elif (thetaPDiff == 0) and (thetaDDiff < 0):
                gamma[indexInstant, indexStride] = -90
            elif (thetaPDiff < 0) and (thetaDDiff == 0):
                gamma[indexInstant, indexStride] = -180
            elif (thetaPDiff == 0) and (thetaDDiff == 0):
                gamma[indexInstant, indexStride] = 0
            if gamma[indexInstant, indexStride] < 0:
                gamma[indexInstant, indexStride] = gamma[indexInstant, indexStride]+360
    gammaMean = np.zeros((1, thetaP.shape[0]-1))
    r = np.zeros((1, thetaP.shape[0] - 1))
    cav = np.zeros((1, thetaP.shape[0] - 1))
    pattern = {
        'in_phase_proximalD_+proximal_+distal': 0,
        'in_phase_distalD_+proximal_+distal': 0,
        'anti_phase_distalD_-proximal_+distal': 0,
        'anti_phase_proximalD_-proximal_+distal': 0,
        'in_phase_proximalD_-proximal_-distal': 0,
        'in_phase_distalD_-proximal_-distal': 0,
        'anti_phase_distalD_+proximal_-distal': 0,
        'anti_phase_proximalD_+proximal_-distal': 0
    }
    for indexInstant in range(thetaP.shape[0]-1):
        xb = np.mean(np.cos(np.deg2rad(gamma[indexInstant, :])))
        yb = np.mean(np.sin(np.deg2rad(gamma[indexInstant, :])))
        if (xb > 0) and (yb > 0):
            gammaMean[0, indexInstant] = np.degrees(math.atan(yb/xb))
        elif xb < 0:
            gammaMean[0, indexInstant] = np.degrees(math.atan(yb/xb))+180
        elif (xb > 0) and (yb < 0):
            gammaMean[0, indexInstant] = np.degrees(math.atan(yb/xb))+360
        elif (xb == 0) and (yb > 0):
            gammaMean[0, indexInstant] = 90
        elif (xb == 0) and (yb < 0):
            gammaMean[0, indexInstant] = -90
        elif (xb == 0) and (yb == 0):
            gammaMean[0, indexInstant] = np.nan
        else:
            gammaMean[0, indexInstant] = 0
        r[0, indexInstant] = np.sqrt(np.power(xb, 2)+np.power(yb, 2))
        cav[0, indexInstant] = np.sqrt((2*(1-r[0, indexInstant])))*(180/np.pi)
        if (gammaMean[0, indexInstant] >= 22.5) and (gammaMean[0, indexInstant] < 67.5):
            pattern['in_phase_proximalD_+proximal_+distal'] = pattern['in_phase_proximalD_+proximal_+distal'] + 1
        elif (gammaMean[0, indexInstant] >= 67.5) and (gammaMean[0, indexInstant] < 112.5):
            pattern['in_phase_distalD_+proximal_+distal'] = pattern['in_phase_distalD_+proximal_+distal'] + 1
        elif (gammaMean[0, indexInstant] >= 112.5) and (gammaMean[0, indexInstant] < 157.5):
            pattern['anti_phase_distalD_-proximal_+distal'] = pattern['anti_phase_distalD_-proximal_+distal'] + 1
        elif (gammaMean[0, indexInstant] >= 157.5) and (gammaMean[0, indexInstant] < 202.5):
            pattern['anti_phase_proximalD_-proximal_+distal'] = pattern['anti_phase_proximalD_-proximal_+distal'] + 1
        elif (gammaMean[0, indexInstant] >= 202.5) and (gammaMean[0, indexInstant] < 247.5):
            pattern['in_phase_proximalD_-proximal_-distal'] = pattern['in_phase_proximalD_-proximal_-distal'] + 1
        elif (gammaMean[0, indexInstant] >= 247.5) and (gammaMean[0, indexInstant] < 292.5):
            pattern['in_phase_distalD_-proximal_-distal'] = pattern['in_phase_distalD_-proximal_-distal'] + 1
        elif (gammaMean[0, indexInstant] >= 292.5) and (gammaMean[0, indexInstant] < 337.5):
            pattern['anti_phase_distalD_+proximal_-distal'] = pattern['anti_phase_distalD_+proximal_-distal'] + 1
        elif (gammaMean[0, indexInstant] >= 337.5) and (gammaMean[0, indexInstant] < 22.5):
            pattern['anti_phase_proximalD_+proximal_-distal'] = pattern['anti_phase_proximalD_+proximal_-distal'] + 1
    return gammaMean, cav, pattern



