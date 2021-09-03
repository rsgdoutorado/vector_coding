import numpy as np
from functions import vector_coding


if __name__ == '__main__':
    x = np.loadtxt('x.txt')
    y = np.loadtxt('y.txt')
    gamma, cav, pattern = vector_coding(np.transpose(x), np.transpose(y))
    for i in pattern.keys():
        print(f'{i}: {pattern.get(i)}')
    print('End...')