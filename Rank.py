import numpy as np
import pandas as pd
import networkx as nx

def geneRank(G, E, s = .85, maxerr = .00001, max_iter=100, weight='weight'):
    """
    Computes the pagerank for each of the n states
    Parameters
    ----------
    G: matrix representing state transitions
       Gij is a binary value representing a transition from state i to j.
    s: probability of following a transition. 1-s probability of teleporting
       to another state.
    maxerr: if the sum of pageranks between iterations is bellow this we will
            have converged.
    """
    if not G.is_directed():
        D = G.to_directed()
    else:
        D = G

    W = nx.stochastic_graph(D, weight=weight)
    N = W.number_of_nodes()
 
    x = dict.fromkeys(W, 1.0 / N)

    p = dict.fromkeys(W, 1.0 / N)
 
    dangling_weights = p

    dangling_nodes = [n for n in W if W.out_degree(n, weight=weight) == 0.0]
 
    for _ in range(max_iter):
        xlast = x
        x = dict.fromkeys(xlast.keys(), 0)
        danglesum = s * sum(xlast[n] for n in dangling_nodes)
        for n in x:
 
            for nbr in W[n]:
                x[nbr] += s * xlast[n] * W[n][nbr][weight]
            x[n] += danglesum * dangling_weights[n] + (1.0 - s) * np.mean(expD.loc[n,])
 
        err = sum([abs(x[n] - xlast[n]) for n in x])
        if err < N*maxerr:
            return x
    raise NetworkXError('generank: power iteration failed to converge '
                        'in %d iterations.' % max_iter)

pinG = pd.read_csv("RPIN.csv")
pinG['weight'] = 1
G = nx.from_pandas_dataframe(pinG,'Source_A','Source_B',['weight'])
expD = pd.read_csv("RGEP.csv")
expD = expD.set_index('Gene')
expD = (expD - expD.min()) / (expD.max() - expD.min())
df = pd.DataFrame.from_dict(geneRank(G,expD),orient='index')
df.to_csv('RankedGenes.csv')
