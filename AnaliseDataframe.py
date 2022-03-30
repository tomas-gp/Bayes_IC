from re import A
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d


# dado = pd.read_csv('csvs/MaiorIntervalo.csv', delimiter=',')

# dado = dado.query('Mean_Abs < 15  &  Std_Corr > 0.85')
# #
# print('Quantidade de conjuntos dado limite: {0}\n'.format(len(dado.index)))
# #
# # print('Min: {0} | Var: {1} | Max: {2}'.format(dado['muP'].min(), dado['muP'].var(), dado['muP'].max()))
# # print('Min: {0} | Var: {1} | Max: {2}'.format(dado['sP'].min(), dado['sP'].var(), dado['sP'].max()))
# # print('Min: {0} | Var: {1} | Max: {2}\n '.format(dado['Nsig'].min(), dado['Nsig'].var(), dado['Nsig'].max()))


# print(dado.sort_values(by=['Mean_Abs'], ascending=[True]))

# fig = plt.figure()
# ax = plt.axes(projection='3d')

# ax.scatter3D(dado['muP'], dado['sP'], dado['Nsig'], c=dado['Mean_Abs'])

# ax.set_xlabel('muP')
# ax.set_ylabel('sP')
# ax.set_zlabel('Nsig')



# # plt.plot(dado['Mean_Abs'], dado['Std_Corr'], '.')

# plt.show()


a = [2, 4, 5]
b = [1, 3, 2]

c = np.array(a) * b
print(c)