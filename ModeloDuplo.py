import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import pandas as pd

'''
Dados:
'''

freq = '400'
bb = 0

azi_cipic = pd.read_csv('csvs/cipic_AzimItd_NYCD1.csv')
fi_cipic = pd.read_csv('csvs/cipic_ItdFi_NYCD1.csv')
beh = pd.read_csv('csvs/behNYC.csv')

itd_beh = beh['itd'].to_numpy()
azimmean_beh = beh[freq + 'm'].to_numpy()
azimstd_beh = beh[freq + 's'].to_numpy()

itdazim_cipic = azi_cipic[freq].to_numpy()
azim_cipic = np.linspace(-90, 90, len(itdazim_cipic))

temp = fi_cipic[freq].isna().to_numpy()
temp = np.concatenate(([0], np.where(temp[1:] ^ temp[:-1])[0] + 1, [temp.size]))

fiitd_cipic = fi_cipic[freq].dropna().to_numpy()
itd_cipic = np.linspace(-940 + temp[1] * 5, -940 + temp[2] * 5 - 5, len(fiitd_cipic))

if bb == 1:
    azimmean_beh1 = beh['600' + 'm'].to_numpy()
    azimstd_beh1 = beh['600' + 's'].to_numpy()

    itdazim_cipic1 = azi_cipic['600'].to_numpy()
    azim_cipic1 = np.linspace(-90, 90, len(itdazim_cipic1))

    temp = fi_cipic['600'].isna().to_numpy()
    temp = np.concatenate(([0], np.where(temp[1:] ^ temp[:-1])[0] + 1, [temp.size]))

    fiitd_cipic1 = fi_cipic['600'].dropna().to_numpy()
    itd_cipic1 = np.linspace(-940 + temp[1] * 5, -940 + temp[2] * 5 - 5, len(fiitd_cipic1))


'''
Funções
'''

def gauss(x, mu, sigma):
    return np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)


def lorentzian(x, x0, gam):
    return gam ** 2 / (gam ** 2 + (x - x0) ** 2)


def pdf_mean_std(xval, pdf):
    X = xval
    W = pdf
    N = len(X)
    Wmean = np.sum(X * W) / np.sum(W)
    Wstd = np.sqrt((np.sum(W * (X - Wmean) ** 2)) / ((N - 1) * np.sum(W) / N))
    return Wmean, Wstd


'''
Visualização dos parametros selecionados
'''

Mp, sP, N = [64, 10, 5]

x = np.linspace(-180, 180, 361)

cs_azim = CubicSpline(itdazim_cipic, azim_cipic)
cs_fi = CubicSpline(itd_cipic, fiitd_cipic)

if bb == 1:
    cs_azim1 = CubicSpline(itdazim_cipic1, azim_cipic1)
    cs_fi1 = CubicSpline(itd_cipic1, fiitd_cipic1)

Std_pos = np.array([])
Std_beh = np.array([])
MeanArr = np.array([])

priorE = lorentzian(x, -Mp, sP)
priorD = lorentzian(x, +Mp, sP)
priorE = priorE / np.trapz(priorE)
priorD = priorD / np.trapz(priorD)

'''
Plot
'''

subplot = np.arange(len(itd_beh))
count = 0
Mean_Abs = 0
for itd in itd_beh:
    behav_gauss = gauss(x, azimmean_beh[itd_beh == itd], azimstd_beh[itd_beh == itd])
    behav_gauss = behav_gauss / np.trapz(behav_gauss)

    likelihood = lorentzian(x, cs_azim(itd), N / cs_fi(itd))
    if bb == 1:
        likelihood1 = lorentzian(x, cs_azim(itd), N / cs_fi1(itd))
        likelihood = likelihood * likelihood1
    likelihood = likelihood / np.trapz(likelihood)

    arr = cs_azim(itd)
    bol = x == round(arr.tolist())
    pos = np.argwhere(bol)

    posterior = priorE[pos] * priorE * likelihood + priorD[pos] * priorD * likelihood
    posterior = posterior / np.trapz(posterior)

    p_mean, p_std = pdf_mean_std(x, posterior)

    MeanArr = np.append(MeanArr, p_mean)

    if bb == 1:
        Mean_Abs = np.abs(p_mean - azimmean_beh1[itd_beh == itd]) + Mean_Abs
    else:
        Mean_Abs = np.abs(p_mean - azimmean_beh[itd_beh == itd]) + Mean_Abs

    if bb == 1:
        Std_pos = np.append(Std_pos, p_std)
        Std_beh = np.append(Std_beh, azimstd_beh1[itd_beh == itd])
    else:
        Std_pos = np.append(Std_pos, p_std)
        Std_beh = np.append(Std_beh, azimstd_beh[itd_beh == itd])

    plt.subplot(3, 4, count + 1)
    plt.title(
        f'Stimulus: {np.round(cs_azim(itd), 2)} Deg \nMean Response: {str(int(azimmean_beh[itd_beh == itd]))} '
        f'Deg ; Model Estimation: {str(round(p_mean, 2))} Deg',
        fontsize=8)
    plt.plot(x, priorE, 'k')
    plt.plot(x, priorD, 'k')
    plt.plot(x, likelihood, 'r')
    plt.plot(x, posterior[0], 'b')
    plt.axvline(p_mean, linestyle='--')
    plt.plot(x, behav_gauss, 'c')
    count = count + 1



print(Mean_Abs)
Std_Corr = np.corrcoef(Std_pos, Std_beh)
print(Std_Corr[0][1])

x = np.linspace(-90, 90, 11)
plt.subplot(3, 4, count + 1)
plt.title('STD do modelo (superior) \n STD do comportamento(inferior)', fontsize=8)
plt.plot(x, Std_beh)
plt.plot(x, Std_pos)
plt.show()
