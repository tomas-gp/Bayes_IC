import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt



'''
Dados:
'''

freq = '400'

# Abrindo os Csvs
azi_cipic = pd.read_csv('csvs/cipic_AzimItd_NYCD1.csv')
# fi_cipic = pd.read_csv('csvs/cipic_ItdFi_NYCD1.csv')
beh = pd.read_csv('csvs/behNYC.csv')
sdCIPIC = pd.read_csv('csvs\sdCIPICcabecamediaNY.csv')

# Criando variaveis do Behaviour
itd_beh = beh['itd'].to_numpy()
azimmean_beh = beh[freq + 'm'].to_numpy()
azimstd_beh = beh[freq + 's'].to_numpy()

# Criando variaveis do ITD -> Azim Cipic
itdazim_cipic = azi_cipic[freq].to_numpy()

itdazim_cipic = -itdazim_cipic[::-1]/2 + itdazim_cipic/2 #


azim_cipic = np.linspace(-90, 90, len(itdazim_cipic))

# Criando variaveis de Fisher (outdated)
# temp = fi_cipic[freq].isna().to_numpy()
# temp = np.concatenate(([0], np.where(temp[1:] ^ temp[:-1])[0] + 1, [temp.size]))
# fiitd_cipic = fi_cipic[freq].dropna().to_numpy()

# Criando variaveis do sd CIPIC (AZIM ou ITD(não implementado))
temp2 = sdCIPIC['sdITD' + freq ].isna().to_numpy()
temp2 = np.concatenate(([0], np.where(temp2[1:] ^ temp2[:-1])[0] + 1, [temp2.size]))
sdAZIM_cipic = sdCIPIC['sdITD' + freq ].dropna().to_numpy()
sdAZIM_cipic = sdAZIM_cipic[::-1]/2 + sdAZIM_cipic/2 #

#sdITD_cipic = sdCIPIC['2sdITD' + freq ].to_numpy() (parte ITD)


# Vetor para usar pro FI/ITD
# itd_cipic = np.linspace(-940 + temp[1] * 5, -940 + temp[2] * 5 - 5, len(fiitd_cipic))


'''
Funções
'''


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
Parâmetros
'''

Mp, sP = [64, 10]

x = np.linspace(-180, 180, 361)

cs_azim = CubicSpline(itdazim_cipic, azim_cipic)
# cs_fi = CubicSpline(itd_cipic, fiitd_cipic)
cs_sdcipic = CubicSpline(azim_cipic, sdAZIM_cipic)


'''
Seek-Until
'''

dataset = np.zeros((181, 1))
# col = ['Prior On', 'Prior Central', 'Prior off']
seek = pd.DataFrame(dataset)#, columns=col)

priorE = lorentzian(x, -Mp, sP)
priorD = lorentzian(x, +Mp, sP)
priorUpdate = 0*priorE
priorE = priorE / np.trapz(priorE)
priorD = priorD / np.trapz(priorD)


count = 0
for itd in itdazim_cipic:
    novo = cs_azim(itd)
    print(count)
    trials = 0
    error = 10
    while True:
        if trials >= 1:
            priorUpdate = posterior

        sd = cs_sdcipic(novo)

        likelihood = lorentzian(x, novo, sd)
        likelihood = likelihood / np.trapz(likelihood) #[float(i) / sum(likelihood) for i in likelihood]

        # Encontra o indice do valor de azimute EX: (-90) = [0]

        pos = np.argwhere(x == round(novo.tolist()))

        # posterior = priorE[pos] * priorE * likelihood + priorD[pos] * priorD * likelihood 
        posterior = priorE[pos] * priorE * likelihood + priorD[pos] * priorD * likelihood + priorUpdate[pos] * priorUpdate* likelihood


        posterior = np.squeeze(np.asarray(posterior))  # Acho que pode excluir
        posterior = posterior / np.trapz(posterior) #[float(i) / sum(posterior) for i in posterior]
        posterior = np.array(posterior) # Acho que pode excluir

        #rand = np.random.choice(x, p=np.squeeze(np.asarray(posterior)))  #Pega valor aleatorio na PDF

        p_mean, p_std = pdf_mean_std(x, posterior)



        # plt.plot(x, priorE, 'k')
        # plt.plot(x, priorD, 'k', label='Prior')
        # plt.plot(x, priorUpdate, 'k', linestyle='--', label='PriorUpdate')

        # plt.plot(x, likelihood, 'r', label='Likelihood')
        # plt.plot(x, posterior, 'b', label='Posterior')
        # plt.axvline(p_mean, linestyle='--', label='Model estimation')
        # plt.axvline(novo, color='r', linestyle='--', label='Estim')
        # plt.xticks(np.arange(-180, 180, step=15))

        # plt.legend()
        # plt.show()

        novo = novo - p_mean

        posterior = np.roll(posterior, int(-p_mean))




        if novo >90:
            print('opa')    
            novo = np.float64(90)
        if novo <-90:
            print('opa')   
            novo = np.float64(-90)


        error = novo
        

        trials = trials + 1

        if np.abs(error) < 3:
            break

    seek.iat[count, 0] = trials
    count = count + 1
# Salvando
seek.to_csv('csvs/seek_repetition_rand3.csv')





# import numpy as np
# import pandas as pd
# from scipy.interpolate import CubicSpline
# import matplotlib.pyplot as plt



# '''
# Dados:
# '''

# freq = '400'

# # Abrindo os Csvs
# azi_cipic = pd.read_csv('csvs/cipic_AzimItd_NYCD1.csv')
# # fi_cipic = pd.read_csv('csvs/cipic_ItdFi_NYCD1.csv')
# beh = pd.read_csv('csvs/behNYC.csv')
# sdCIPIC = pd.read_csv('csvs\sdCIPICcabecamediaNY.csv')

# # Criando variaveis do Behaviour
# itd_beh = beh['itd'].to_numpy()
# azimmean_beh = beh[freq + 'm'].to_numpy()
# azimstd_beh = beh[freq + 's'].to_numpy()

# # Criando variaveis do ITD -> Azim Cipic
# itdazim_cipic = azi_cipic[freq].to_numpy()
# azim_cipic = np.linspace(-90, 90, len(itdazim_cipic))

# # Criando variaveis de Fisher (outdated)
# # temp = fi_cipic[freq].isna().to_numpy()
# # temp = np.concatenate(([0], np.where(temp[1:] ^ temp[:-1])[0] + 1, [temp.size]))
# # fiitd_cipic = fi_cipic[freq].dropna().to_numpy()

# # Criando variaveis do sd CIPIC (AZIM ou ITD(não implementado))
# temp2 = sdCIPIC['sdITD' + freq ].isna().to_numpy()
# temp2 = np.concatenate(([0], np.where(temp2[1:] ^ temp2[:-1])[0] + 1, [temp2.size]))
# sdAZIM_cipic = sdCIPIC['sdITD' + freq ].dropna().to_numpy() 
# #sdITD_cipic = sdCIPIC['2sdITD' + freq ].to_numpy() (parte ITD)


# # Vetor para usar pro FI/ITD
# # itd_cipic = np.linspace(-940 + temp[1] * 5, -940 + temp[2] * 5 - 5, len(fiitd_cipic))


# '''
# Funções
# '''


# def lorentzian(x, x0, gam):
#     return gam ** 2 / (gam ** 2 + (x - x0) ** 2)


# def pdf_mean_std(xval, pdf):
#     X = xval
#     W = pdf
#     N = len(X)
#     Wmean = np.sum(X * W) / np.sum(W)
#     Wstd = np.sqrt((np.sum(W * (X - Wmean) ** 2)) / ((N - 1) * np.sum(W) / N))

#     return Wmean, Wstd


# '''
# Parâmetros
# '''

# Mp, sP = [64, 10]

# x = np.linspace(-180, 180, 361)

# cs_azim = CubicSpline(itdazim_cipic, azim_cipic)
# # cs_fi = CubicSpline(itd_cipic, fiitd_cipic)
# cs_sdcipic = CubicSpline(azim_cipic, sdAZIM_cipic)


# '''
# Seek-Until
# '''
# repetitionsss = 1
# dataset = np.zeros((181, repetitionsss))
# # col = ['Prior On', 'Prior Central', 'Prior off']
# seek = pd.DataFrame(dataset)#, columns=col)

# priorE = lorentzian(x, -Mp, sP)
# priorD = lorentzian(x, +Mp, sP)
# priorUpdate = 0*priorE
# priorE = priorE / np.trapz(priorE)
# priorD = priorD / np.trapz(priorD)



# for repetition in range(repetitionsss):
#     count = 0
#     for itd in itdazim_cipic:
#         novo = cs_azim(itd)
#         print(count)
#         trials = 0
#         error = 10
#         while True:
#             if trials >= 1:
#                 priorUpdate = posterior

#             sd = cs_sdcipic(novo)

#             likelihood = lorentzian(x, novo, sd)
#             likelihood = likelihood / np.trapz(likelihood) #[float(i) / sum(likelihood) for i in likelihood]

#             # Encontra o indice do valor de azimute EX: (-90) = [0]
#             bol = x == round(novo.tolist())
#             pos = np.argwhere(bol)

#             # posterior = priorE[pos] * priorE * likelihood + priorD[pos] * priorD * likelihood 
#             posterior = priorE[pos] * priorE * likelihood + priorD[pos] * priorD * likelihood + priorUpdate[pos] * priorUpdate* likelihood


#             posterior = np.squeeze(np.asarray(posterior))
#             posterior = posterior / np.trapz(posterior) #[float(i) / sum(posterior) for i in posterior]
#             posterior = np.array(posterior)

#             #rand = np.random.choice(x, p=np.squeeze(np.asarray(posterior)))  #Pega valor aleatorio na PDF

#             p_mean, p_std = pdf_mean_std(x, posterior)



#             # plt.plot(x, priorE, 'k')
#             # plt.plot(x, priorD, 'k', label='Prior')
#             # plt.plot(x, priorUpdate, 'k', linestyle='--', label='PriorUpdate')

#             # plt.plot(x, likelihood, 'r', label='Likelihood')
#             # plt.plot(x, posterior, 'b', label='Posterior')
#             # plt.axvline(p_mean, linestyle='--', label='Model estimation')
#             # plt.axvline(novo, color='r', linestyle='--', label='Estim')
#             # plt.xticks(np.arange(-180, 180, step=15))

#             # plt.legend()
#             # plt.show()

#             novo = novo - p_mean

#             posterior = np.roll(posterior, int(-p_mean))




#             if novo >90:
#                 print('opa')    
#                 novo = np.float64(90)
#             if novo <-90:
#                 print('opa')   
#                 novo = np.float64(-90)


#             error = novo
            

#             trials = trials + 1

#             if np.abs(error) < 3:
#                 break

#         seek.iat[count, repetition] = trials
#         count = count + 1
# # Salvando
# seek.to_csv('seek_repetition_rand.csv')
