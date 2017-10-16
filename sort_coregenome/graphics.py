#!/usr/bin/python
def graph(df,cutoff):
    import matplotlib.pyplot as plt
    import numpy as np
    #rng = np.random.RandomState(10)  # deterministic random data


    fig, axes = plt.subplots(nrows=3)
    fig.subplots_adjust(hspace=0.4)
    ax0, ax1, ax2= axes.flatten()

    
    ax0.hist(df["hitsIn"], bins=35,color="darkgreen")  # arguments are passed to np.histogram
    ax0.set_title("Histogram of gene distribution of Out-group")
   
    
    ax1.hist(df["hitsOut"], bins=10,color="teal")  # arguments are passed to np.histogram
    ax1.set_title("Histogram of gene distribution of out-group")
    

    hb = ax2.hexbin(x=df['hitsIn'], y=df['hitsOut'],bins='log',gridsize=30,cmap="YlGn")
    cb = fig.colorbar(hb)
    cb.set_label('log10(N)')
    plt.show()