import matplotlib.pyplot as plt

import math
import random
import numpy as np
import pandas as pd

n_samples_low = 100000


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

    b = 0
    for i in range(samples):
        b = random.random()*(i1-i0) + i0
        b = func(s)*(i1-i0)
        b = (b*(i) + s)/(i+1)


    def mutate(x):
        x = x + (random.random()-0.5)*(i1-i0)/10.0;
        if(x>i1):
            x = i1-x+i1
        if(x<i0):
            x = i0-x+i0

        return [x, 10.0/(i1-i0)]

    def reset():
        return [random.random()*(i1-i0) + i0, 1.0/(i1-i0)];


    [x,px] = reset()
    v = func(x)
    first_v = v
    first_p = px
    samples = n_samples_low
    sum = 0
    plot_data = [0]*samples
    sample_pts = [0]*samples
    for i in range(samples):

        ecs = random.random()
        if(ecs>0.01):
            [xi,pxi] = reset()
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
        a = v
        sum = (sum*i + a)/(i+1)
        plot_data[i] = sum/b
        sample_pts[i] = x
    return plot_data,sample_pts




def guidedMC(func, i0,i1):
    pass

def makeHistogram(i0,i1,n,samples):
    count = [0]*n
    xis = np.arange(i0,i1-0.0001,(i1-i0)/n)

    for s in samples:
        count[int(s*n - 0.0001)] +=1/len(samples)

    return [xis,count]




def main():
    i0 = 0
    i1 = 1
    x = np.arange(0,n_samples_low,1)
    a = montecarlo(f,i0,i1)
    a = np.array(a)
    fig = plt.figure()

    fig.add_subplot(2, 3, 1)
    plt.plot(x,a)


    [a,samples] = mcmc(f,i0,i1,b)

    fig.add_subplot(2, 3, 2)
    plt.plot(x,a)
   # plt.plot(x,a)
    #guidedMC(f,0,1)

    histogram = makeHistogram(i0,i1,50,samples)

    fig.add_subplot(2, 3, 3)
    plt.plot(histogram[0],histogram[1])
    #pd.DataFrame(samples).plot.density()


    xs = np.arange(0,1,0.001)
    ys = np.array(list(map(f,xs)))


    fig.add_subplot(2, 3, 4)
    plt.plot(xs,ys)
    plt.show()


main();