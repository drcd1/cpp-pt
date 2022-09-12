import matplotlib.pyplot as plt

import math
import random
import numpy as np
import pandas as pd

n_samples_low = 1000


def gauss(x, m, std):
    return (1.0/(std*math.sqrt(2*math.pi) ))*math.exp(-(x-m)*(x-m)*0.5/(std*std))


def f(x):
    return gauss(x,0,0.1) + gauss(x,0.6,0.3) + gauss(x,0.8,0.1) + gauss(x,0.2,0.05) + gauss(x,0.5,0.03)

def montecarlo(func, i0,i1):
    samples = n_samples_low
    sum = 0
    plot_data = [0]*samples
    for i in range(samples):
        s = random.random()*(i1-i0) + i0
        s = func(s)*(i1-i0)
        sum = (sum*(i) + s)/(i+1)
        plot_data[i] = sum

    return plot_data





def mcmc(func,i0,i1):
    def mutate(x):
        x = x + (random.random()-0.5)*(i1-i0)/10.0;
        return [x, 10.0/(i1-i0)]

    def reset():
        return [random.random()*(i1-i0) + i0, 1.0/(i1-i0)];

    [x,px] = reset()
    v = func(x)
    samples = n_samples_low
    sum = 0
    plot_data = [0]*samples
    sample_pts = [0]*samples
    for i in range(samples):

        if(i%50==0):
            [x,px] = reset()
            v = func(x)
        else:
            [xi,pxi] = mutate(x)

            vi = func(xi)
            alpha = vi/v

            if(x<i0 or x>i1):
                alpha = 0

            if(alpha>1.0):
                px = px*pxi
                x = xi
                v = vi
            else:
                ecs = random.random();
                if(ecs<alpha):
                    px = px*pxi*alpha
                    x = xi
                    v = vi
                else:
                    px = px*pxi*(1.0-alpha)
        a = v/px;
        sum = (sum*i + a)/(i+1)
        plot_data[i] = sum
        sample_pts[i] = x
    return plot_data,sample_pts




def guidedMC(func, i0,i1):
    pass





def main():
    x = np.arange(0,n_samples_low,1)
    a = montecarlo(f,0,1)
    a = np.array(a)
    fig = plt.figure()

    fig.add_subplot(2, 3, 1)
    plt.plot(x,a)

    [a,samples] = mcmc(f,0,1)

    fig.add_subplot(2, 3, 2)
    plt.plot(x,a)
   # plt.plot(x,a)
    #guidedMC(f,0,1)


    fig.add_subplot(2, 3, 3)
    pd.DataFrame(samples).plot.density()


    xs = np.arange(0,1,0.001)
    ys = np.array(list(map(f,xs)))


    fig.add_subplot(2, 3, 4)
    plt.plot(xs,ys)
    plt.show()


main();