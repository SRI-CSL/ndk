from __future__ import print_function

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
# %matplotlib inline
import matplotlib.pyplot as plt

import seaborn as sns


def pca_reduce_waveforms(waveforms, ncomp=3):
    col_labels = ['samp'+str(k) for k in range(len(waveforms[0]))]
    df = pd.DataFrame(waveforms, columns=col_labels)
    pca = PCA(n_components=ncomp)
    pca_result = pca.fit_transform(df[col_labels].values)
    return pca_result


def pca_plot_waveforms(waveforms, ncomp=3, use_tsne=False):
    col_labels = ['samp'+str(k) for k in range(len(waveforms[0]))]
    df = pd.DataFrame(waveforms, columns=col_labels)
    pca = PCA(n_components=ncomp)
    pca_result = pca.fit_transform(df[col_labels].values)
    df['pca-one']   = pca_result[:,0]
    df['pca-two']   = pca_result[:,1]
    df['pca-three'] = pca_result[:,2]
    df['y'] = [k for k in range(len(waveforms))]
    df['label'] = df['y'].apply(lambda i: str(i))
    
    np.random.seed(42)
    rndperm = np.random.permutation(df.shape[0])

    print('Explained variation for each principal component: {}'.format(pca.explained_variance_ratio_))

    plt.figure(figsize=(16,10))
    sns.scatterplot(x="pca-one", 
                    y="pca-two", 
                    hue="pca-three", 
                    #palette=sns.color_palette("hls", 10),
                    data=df.loc[rndperm,:],
                    legend="full",
                    alpha=0.3)

    if use_tsne:
        tsne = TSNE(n_components=2, verbose=1, perplexity=30, n_iter=500)
        tsne_results = tsne.fit_transform(df[col_labels].values)

        df['tsne-one'] = tsne_results[:,0]
        df['tsne-two'] = tsne_results[:,1]

        plt.figure(figsize=(16,10))
        sns.scatterplot(x="tsne-one", y="tsne-two", data=df.loc[rndperm,:], legend="full", alpha=0.3)

    return df
