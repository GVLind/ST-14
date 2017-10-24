#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def graph(df, cutoff):
    
    #rng = np.random.RandomState(10)  # deterministic random data


    fig, axes = plt.subplots(nrows = 2)
    fig.subplots_adjust(hspace = 0.4, right = 0.87,left = 0.18)
    ax0, ax1 = axes.flatten()

    
    ax0.hist(df["hitsIn"], bins = 75, color = "black")  # arguments are passed to np.histogram
    ax0.set_title("Histogram of gene distribution of In-group")
   
    
    ax1.hist(df["hitsOut"], bins = 10, color = "black")  # arguments are passed to np.histogram
    ax1.set_title("Histogram of gene distribution of out-group")
    
    plt.show()

    plt.hexbin(x = df['hitsIn'], y = df['hitsOut'], bins = "log", gridsize = 30, cmap = "CMRmap_r")
    cb = plt.colorbar()
    cb.set_label('No strains log10(N)')
    plt.show()


if __name__ == "__main__":
    cutoff = [0.95, 0.5]
    df = pd.read_csv("sample.csv", sep = "\t")
    graph(df, cutoff)
    #cmap = "BuPu, RdBu_r, binary, cubehelix_r, CMRmap_r"