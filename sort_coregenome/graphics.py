#!/usr/bin/python
def graph(df,cutoff):
    import matplotlib.pyplot as plt
    import numpy as np
    #rng = np.random.RandomState(10)  # deterministic random data


    fig, axes = plt.subplots(nrows=2)
    fig.subplots_adjust(hspace=0.4)
    ax0, ax1= axes.flatten()

    
    ax0.hist(df["hitsIn"], bins=35,color="darkgreen")  # arguments are passed to np.histogram
    ax0.set_title("Histogram of gene distribution of In-group")
   
    
    ax1.hist(df["hitsOut"], bins=10,color="teal")  # arguments are passed to np.histogram
    ax1.set_title("Histogram of gene distribution of out-group")
    
    plt.show()

    plt.hexbin(x=df['hitsIn'], y=df['hitsOut'],bins='log',gridsize=50,cmap="Reds")
    cb = plt.colorbar()
    cb.set_label('log10(N)')
    plt.show()