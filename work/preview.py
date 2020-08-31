from os.path import isfile
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy import units as u
import lightkurve as lk

# PATH = "Data/TESS/GI 11047/r1/spoc/"
# CATALOGUE="TIC"
# df = pd.read_excel("Correspondence/TESS_results_excel.xlsx", skiprows=1)
# df = df.rename(columns={df.columns[0]: 'ID', df.columns[-1]: 'flag'})

# m = np.array([isinstance(i, str) and i[0] in ['Y', 'y'] for i in df['flag']])
# Δν_list = df[m][['ID', 'numax', 'numax_error', 'dnu', 'dnu_error']].rename(columns={'numax': 'ν_max', 'dnu': 'Δν'})

# files = [f"{PATH}/{int(_)}_psd_tot_nw.txt" for _ in Δν_list['ID']]

# def job(r):
#     ID = int(r['ID'])
#     # file = f"{PATH}/{ID}_psd_tot_nw.txt"  
#     Δν = float(r['Δν'])
#     ν_max = float(r['ν_max'])

# PATH = "Kepler"
# CATALOGUE="KIC"
# files = sorted(glob(f"{PATH}/*-psd_slc*.pow"))
# Δν_list = pd.read_csv(f"{PATH}/Δν", delim_whitespace=False, skiprows=1, names=['ID', 'ν_max', 'e_ν_max', 'Δν', 'e_Δν'])

# from astropy.io import fits
# app = fits.getdata("Seed/Appourchaux 2012.fits")
# app_KIC_all = app.KIC.astype(int)
# app_KIC = np.unique(app_KIC_all)

# PATH = "Data/C17"
# CATALOGUE = "EPIC"
# files = sorted(glob(f"{PATH}/*-psd_slc*.pow"))
# Δν_list = pd.read_csv(f"{PATH}/Δν", delim_whitespace=True, skiprows=1,
#     names=['ID', 'Δν', 'ν_max'])

# def job(file):
#     # ID = int(file.split("kplr")[1].split("_")[0])
#     ID = int(file.split("ktwo")[1].split("_")[0])
#     r = Δν_list[Δν_list['ID'] == ID]
#     Δν = float(r['Δν'])
#     ν_max = float(r['ν_max'])

project = "dublicates"

if project == "TESS C17":
    PATH = "Data/TESS/GI 11047/r1/spoc/"
    CATALOGUE="TIC"
    df = pd.read_excel("Correspondence/TESS_results_excel.xlsx", skiprows=1)
    df = df.rename(columns={df.columns[0]: 'ID', df.columns[-1]: 'flag'})

    m = np.array([isinstance(i, str) and i[0] in ['Y', 'y'] for i in df['flag']])
    Δν_list = df[m][['ID', 'numax', 'numax_error', 'dnu', 'dnu_error']].rename(columns={'numax': 'ν_max', 'dnu': 'Δν'})
    Δν_list.to_csv("TESS_Δν.csv", index=False)

    files = [f"{PATH}/{int(_)}_psd_tot_nw.txt" for _ in Δν_list['ID']]

elif project == "K2 subgiants":
    PATH = "Data/C17"
    CATALOGUE = "EPIC"
    files = sorted(glob(f"{PATH}/*.pow"))
    Δν_list = pd.read_csv(f"{PATH}/Δν", delim_whitespace=True, skiprows=1, names=['ID', 'Δν', 'ν_max'])

elif project == "dublicates":
    PATH = "Data/Dublicates/"
    CATALOGUE="EPIC"
    files = sorted(glob(f"{PATH}/*/*_psd_tot_nw.txt"))
    Δν_list = pd.read_csv(f"{PATH}/Δν", delim_whitespace=False, skiprows=1, names=['ID', 'Δν', 'e_Δν', 'ν_max', 'e_νmax'])


# def job(file):

#     if CATALOGUE == "KIC":
#         ID = int(file.split("kplr")[1].split("_")[0])
#     elif CATALOGUE == "EPIC":
#         ID = int(file.split("ktwo")[1].split("_")[0])
#     elif CATALOGUE == "TIC":
#         ID = int(basename(file).split('_')[0])
#     r = Δν_list[Δν_list['ID'] == ID]

def job(r, save=True):
    ID = int(r['ID'])
    if CATALOGUE == "TIC":
        file = f"{PATH}/{ID}_psd_tot_nw.txt"
    elif CATALOGUE == "EPIC":
        file = f"{PATH}/{ID}/{ID}_psd_tot.txt"

    Δν = float(r['Δν'])
    ν_max = float(r['ν_max'])
    if Δν < 0:
        return

    # Read in power spectrum
    try:
        file = f"DIAMONDS/data/{CATALOGUE}{ID}.txt"
        ν, ps = np.loadtxt(file).T
    except (FileNotFoundError, OSError):
        # tqdm.write(f"PS for {CATALOGUE} {ID} not found")
        return
    except ValueError:
        print(f"PS for {CATALOGUE} {ID} malformed")
        return

    # Create periodogram object
    pg = lk.Periodogram(ν * u.uHz, u.Quantity(ps), label=f'{CATALOGUE} {ID}')

    # Make echelle diagram
    plt.close('all')
    ss = pg.to_seismology()
    ss.numax = ν_max
    ss.deltanu = Δν
    ss.plot_echelle(smooth_filter_width=Δν/100)

    # get manual mode ID
    if isfile(f'Seed/{CATALOGUE}{ID}.pkl'):
        seed = pd.read_pickle(f'Seed/{CATALOGUE}{ID}.pkl')
        modes = seed.nu_med.values

        plt.scatter(modes % Δν, (modes // Δν) * Δν, c='green')

    # get Appourchaux frequencies
    if CATALOGUE == "KIC" and ID in app_KIC:
        modes = app[app_KIC_all == ID].Freq.astype('float')
        plt.plot(modes % Δν, (modes // Δν) * Δν, 'v', c='red', fillstyle='none')

    # get MLE estimate
    if isfile(f'mle/{CATALOGUE}{ID}.npy'):
        params = np.load(f'mle/{CATALOGUE}{ID}.npy', allow_pickle=True)
        cov = np.array([(lambda J: np.linalg.inv(J.T @ J))(j['jac']) for j in params])
        ν0 = np.array([j['x'][0] for j in params])
        e_ν0 = np.array([np.sqrt(C[0, 0]) for C in cov])

        plt.errorbar(ν0 % Δν, (ν0 // Δν) * Δν, fmt='^',
                        color='orange', fillstyle='none')
    
    # Get DIAMONDS results
    def get_results(n):
        result = f"DIAMONDS/results/{CATALOGUE}{ID}/pb/{n}/peakbagging_parameterSummary.txt"
        if not isfile(result):
            return None

        res = np.loadtxt(result)
        return res.reshape((len(res)//3, 3, -1))

    for i in range(3):
        res = get_results(i)
        if res is None:
            break

        for mode in res:
            ν = mode[0, 1]
            eν1 = mode[0, 4]
            eν2 = mode[0, 5]
            err = [[ν - eν1], [eν2 - ν]]
            plt.errorbar([ν % Δν], [(ν // Δν) * Δν], xerr=err, yerr=err, fmt='.', c='red')
            print(ν)

    if save:
        plt.savefig(f"DIAMONDS/preview/{CATALOGUE}{ID}.png")
    else:
        plt.show()

if __name__ == "__main__":
    list(map(lambda p: job(p[1]), Δν_list.iterrows()))
    # list(map(job, files))
