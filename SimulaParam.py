import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import time

'''
Dados:
'''

freq = '400'

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

itd_cipic = np.linspace(-940 + temp[1]*5, -940 + temp[2]*5 - 5, len(fiitd_cipic))


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
Parâmetros || A definir
'''

# Sigma Priors
sigma_priorE = np.linspace(0, 30, 61)  # 19 - 5:5  ---  37 - 2.5:2.5


# Mean Priors
mu_priorE = np.linspace(0, 80, 161)  # 19 - 5:5  ---  37 - 2.5:2.5


# FreeParameter para o sigma likelihood
num_likelihood = np.linspace(0, 10, 101)  # 51 - 0.2:0.2


'''
Simulação
'''

# Valores fixos
x = np.linspace(-180, 180, 361)
count = 0

cs_azim = CubicSpline(itdazim_cipic, azim_cipic)
cs_fi = CubicSpline(itd_cipic, fiitd_cipic)

lin = len(sigma_priorE)*len(mu_priorE)*len(num_likelihood)
dataset = np.zeros((lin, 5))
col = ['muP', 'sP', 'Nsig', 'Mean_Abs', 'Std_Corr']

data = pd.DataFrame(dataset, columns=col)

start = time.time()  # 1000 trials = 1.45 segs aprox

for muP in mu_priorE:
    print('Mup: ' + str(muP) + ' || Tempo:' + str(round((time.time() - start), 2)))
    for sP in sigma_priorE:
        priorE = lorentzian(x, -muP, sP)
        priorD = lorentzian(x, +muP, sP)
        priorE = priorE / np.trapz(priorE)
        priorD = priorD / np.trapz(priorD)

        for Nsig in num_likelihood:
            Mean_Abs = 0
            Std_pos = np.array([])
            Std_beh = np.array([])
            for itd in itd_beh:
                likelihood = lorentzian(x, cs_azim(itd), Nsig / cs_fi(itd))
                likelihood = likelihood / np.trapz(likelihood)

                arr = cs_azim(itd)
                bol = x == round(arr.tolist())
                pos = np.argwhere(bol)

                posterior = priorE[pos] * priorE * likelihood + priorD[pos] * priorD * likelihood
                posterior = posterior / np.trapz(posterior)

                p_mean, p_std = pdf_mean_std(x, posterior)

                Mean_Abs = np.abs(p_mean - azimmean_beh[itd_beh == itd]) + Mean_Abs
                Std_pos = np.append(Std_pos, p_std)
                Std_beh = np.append(Std_beh, azimstd_beh[itd_beh == itd])
            Std_Corr = np.corrcoef(Std_pos, Std_beh)

            data.iat[count, 0] = muP
            data.iat[count, 1] = sP
            data.iat[count, 2] = Nsig
            data.iat[count, 3] = Mean_Abs
            data.iat[count, 4] = Std_Corr[0][1]

            count = count + 1

data.to_csv('MaiorIntervalo.csv')  # Salva em CSV os dados (Analise em outro código)

