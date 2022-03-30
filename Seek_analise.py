import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

dado = pd.read_csv('csvs\seek_repetition_rand2.csv', delimiter=',')
# dado['mean'] = dado.mean(axis=1)
# print(dado['mean'])
# plt.plot(np.linspace(-90,90,181), dado['mean'], '|')
plt.plot(np.linspace(-90,90,181), dado['0'].tolist(), '-')
plt.show()