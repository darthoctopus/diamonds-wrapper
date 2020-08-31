#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib as mpl
mpl.rcParams['backend'] = 'agg'
import pandas as pd
import matplotlib.pyplot as plt
from os import system
from os.path import isfile, basename

from glob import glob
from tqdm.auto import tqdm


# In[2]:


from astropy.io import fits
app = fits.getdata("Seed/Appourchaux 2012.fits")
app_KIC_all = app.KIC.astype(int)
app_KIC = np.unique(app_KIC_all)


# In[3]:

# PATH = "TESS/"
# CATALOGUE="TIC"
# files = sorted(glob(f"{PATH}/*.txt"))

# PATH = "Kepler"
# CATALOGUE="KIC"
# files = sorted(glob(f"{PATH}/*-psd_slc*.pow"))

project = "TESS C17"

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



#Δν_list = pd.read_csv(f"{PATH}/Δν", delim_whitespace=False, skiprows=1, names=['ID', 'ν_max', 'e_ν_max', 'Δν', 'e_Δν'])
# print(len(files))


# In[5]:


# file = files[0]

# def job(file):

#     if CATALOGUE == "KIC":
#         ID = int(file.split("kplr")[1].split("_")[0])
#     elif CATALOGUE == "EPIC":
#         ID = int(file.split("ktwo")[1].split("_")[0])
#     elif CATALOGUE == "TIC":
#         ID = int(basename(file).split('_')[0])
#     r = Δν_list[Δν_list['ID'] == ID]
#     df = pd.read_csv(file, delim_whitespace=True, skiprows=11, names=['ν','P'])


def job(r):
    ID = int(r['ID'])
    if CATALOGUE == "TIC":
        file = f"{PATH}/{ID}_psd_tot_nw.txt"
    elif CATALOGUE == "EPIC":
        file = f"{PATH}/{ID}/{ID}_psd_tot.txt"
    df = pd.read_csv(file, delim_whitespace=True, names=['ν','P1', 'P2', 'P'])


    Δν = float(r['Δν'])
    νmax = float(r['ν_max'])
    if νmax < 0:
        return


    # plt.plot(df['ν'], df['P'])
    # plt.gca().set_xscale('log')
    # plt.gca().set_yscale('log')
    # plt.axvline(νmax, c='red')

    # def Harvey(ν, t, g):
    #     return np.power(10, t) / (1 + (ν/np.power(10, g))**2)

    def Harvey(ν, A, ν0):
        return 2 * np.sqrt(2) / np.pi * A**2 / (ν0 * (1 + (ν / ν0)**4))

    from scipy.optimize import least_squares

    def model(q):
        p = q[1:]
        return q[0] + np.sum([Harvey(df.ν.values, *p[2*i:2*i+2]) for i in range(len(p)//2)], axis=0)

    def cost(p):
        return (df.P.values - model(p))#/np.sqrt(df.P.values)

    p0 = [0, 10, 100, 10, 500]
    j = least_squares(cost, p0)

    # for i in range(len(p0[1:])//2):
    #     plt.plot(df['ν'].values, Harvey(df['ν'].values, *j['x'][2*i+1:2*i+3]))
    # plt.plot(df['ν'].values, model(j['x']))
    # # plt.plot(df['ν'].values, model(session.bg[:-3]))
    # plt.show()


    # In[ ]:


    from os.path import isfile
    # if isfile(f"./Seed/{CATALOGUE}{ID}.pkl"):
    #     seed = pd.read_pickle(f"Seed/{CATALOGUE}{ID}.pkl")
    #     # for l, r in seed.groupby('ell'):
    #     #     plt.plot(r.nu_med % Δν, r.nu_med, 'o')
        
    #     if CATALOGUE == "KIC" and ID in app_KIC:
    #         modes = app[app_KIC_all == ID]
    #     #     for l in [0, 1, 2, 3]:
    #     #         ml = modes[modes['degree'] == l]
    #     #         ν = ml.Freq.astype(float)
    #             # plt.plot(ν % Δν, ν, '^')
    #     # plt.show()


    # In[7]:


    ε = 1e-10

    w = νmax**0.88 * 0.66

    v = np.concatenate((np.abs(j['x']) * ε, [1, νmax, w]))
    e = np.concatenate((np.ones_like(j['x']) * 3/ε, [100, 1.1, 2]))


    # In[8]:


    j['x']


    # In[10]:


    from diamonds import DiamondsSession

    system(f"mkdir -p 'DIAMONDS/temp/{CATALOGUE}'")
    temp = f"DIAMONDS/temp/{CATALOGUE}/{ID}_PS.txt"
    np.savetxt(temp, np.array([df['ν'].values, df['P'].values]).T)
    session = DiamondsSession(int(ID), temp, n_head=0, CATALOGUE=CATALOGUE)
    # session.bg = 1
    session.bg = v
    session.write_bg_hyperparams(e)
    session.run_bg()


    if isfile(f"Seed/{CATALOGUE}{ID}.pkl"):
        seed = pd.read_pickle(f"Seed/{CATALOGUE}{ID}.pkl")
        modes = seed.nu_med.values

    else:
        return

    if CATALOGUE == "KIC" and ID in app_KIC:
        modes = app[app_KIC_all == ID].Freq.astype('float')

    # WIDTH = 8 # μHz
    WIDTH = Δν/10

    if not isfile(f'mle/{CATALOGUE}{ID}.npy'):
        ν, ps = np.loadtxt(f"DIAMONDS/data/{CATALOGUE}{ID}.txt").T

        from scipy.optimize import least_squares
        def Lorentzian(x, x0, A, Γ):
            # Same parameterisation as DIAMONDS
            return A**2 / ( (x-x0)**2 + (Γ/2)**2 )

        def get_params(ν0):
            m = (ν > ν0 - WIDTH/2) & (ν < ν0 + WIDTH/2)
            j = least_squares(lambda p: (ps[m] - Lorentzian(ν[m], *p)), [ν0, np.max(ps[m]), WIDTH])
            # plt.plot(ν[m], ps[m])
            # plt.plot(ν[m], Lorentzian(ν[m], *j['x']))
            # plt.show()
            return j

        params = np.array([get_params(ν0) for ν0 in modes])
        np.save(f'mle/{CATALOGUE}{ID}', params)

    else:
        params = np.load(f'mle/{CATALOGUE}{ID}.npy', allow_pickle=True)

    cov = np.array([(lambda J: np.linalg.inv(J.T @ J))(j['jac']) for j in params])

    ν0 = np.array([j['x'][0] for j in params])
    e_ν0 = np.array([np.sqrt(C[0, 0]) for C in cov])

    amp = np.abs([j['x'][1] for j in params])
    e_amp = np.abs([np.sqrt(C[1,1]) for C in cov])

    width = np.abs([j['x'][2] for j in params])
    e_width = np.abs([np.sqrt(C[2,2]) for C in cov])

    MIN_SNR = 10
    try:
        m = (e_ν0 < 2.5 * width) & (amp > MIN_SNR * e_amp) & (np.abs(ν0 - modes) < WIDTH)
    except ValueError:
        tqdm.write(f"Array length mismatch for target {CATALOGUE} {ID}")
        return

    # if ID == 212586030:
    #     m_debug = (e_ν0 < width)
    #     print(e_ν0)
    #     print(width)
    #     plt.plot(ν0 % Δν, ν0, 'or')
    #     plt.plot(ν0[m_debug] % Δν, ν0[m_debug], 'o')
    #     plt.show()

    ν0 = ν0[m]
    e_ν0 = e_ν0[m]
    amp = amp[m]
    e_amp = e_amp[m]
    width = width[m]
    e_width = e_width[m]

    l = seed['ell'].values[m]
    
    # plt.plot(ν, ps)
    # plt.scatter(modes, amp)
    # plt.show(block=True)

    # m = amp > 100

    # modes = modes[m]
    # amp = amp[m]

    NMODES = 20
    n_runs = (len(ν0) - 1) // NMODES + 1

    OUTPUT = np.ones((len(ν0), 10)) * np.nan
    OUTPUT[:, 0] = l
    for i in range(n_runs):
        m = np.arange(NMODES) + i * NMODES
        m = m[m < len(ν0)]
        session.pb = ν0[m]
        try:
            session.write_pb_hyperparams(np.maximum(Δν/20, e_ν0[m] * 5), amp[m], width[m],
                e_amp=np.maximum(10, e_amp[m]*20),
                e_width=np.maximum(Δν/20, e_width[m]*5), run=i)
            session.run_pb(run=i)
            result = session._pb_result

            OUTPUT[m, 1] = session._pb_result[::3, 1]
            OUTPUT[m, 2] = session._pb_result[::3, 5]
            OUTPUT[m, 3] = session._pb_result[::3, 4]

            OUTPUT[m, 4] = session._pb_result[1::3, 2]
            OUTPUT[m, 5] = session._pb_result[1::3, 5]
            OUTPUT[m, 6] = session._pb_result[1::3, 4]

            OUTPUT[m, 7] = session._pb_result[2::3, 2]
            OUTPUT[m, 8] = session._pb_result[2::3, 5]
            OUTPUT[m, 9] = session._pb_result[2::3, 4]

        except Exception as e:
            print(amp / e_amp)
            tqdm.write(f"Run {i} failed for {file}: {e}")
    np.savetxt(f"Output/{CATALOGUE}_{ID}.txt", OUTPUT)

from mesa_tricks.utils import tqdm_parallel_map
tqdm_parallel_map(job, [r for _, r in Δν_list.iterrows()], nthreads=5)
# tqdm_parallel_map(job, files, nthreads=6)


