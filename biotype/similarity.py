# Compute similarity between peak arrays, where a peak is in the form <m/z, intensity>

import math
import biotype.peak_extract as pke

def similar_to(pk1, pk2, method='jaccard', rank=5):
    if method == 'rank':
        return rank_similarity(pk1, pk2, good_with=rank)
    if method == 'weighted':
        return rank_similarity(pk1, pk2, good_with=rank, weighted=True)
    else:
        return jaccard_similarity(pk1, pk2)

# def jaccard_similarity(pk1, pk2, d=1):
def jaccard_similarity(pk1, pk2, delta=3):
    """ Compute jaccard similarity of peak arrays of `pk1`, `pk2`
      where two peaks coincidently occurring at the same m/z location are considered similar
    @parms 
      # d: number of digits to round up [NOT USED]
      delta: threshold for idential peaks [NEW]
    """
    # npk1, npk2 = _peak(pk1, g=d), _peak(pk2, g=d)
    npk1, npk2 = pke.align(pk1, pk2, delta)
    common = npk1.keys() & npk2.keys()
    union  = npk1.keys() | npk2.keys()
    
    # n_common: # of co-occurent peaks
    n_common, n_union = len(common), len(union)
    return n_common / float(n_union)

# def rank_similarity(pk1, pk2, d=1, good_with=2, weighted=False):
def rank_similarity(pk1, pk2, delta=3, good_with=2, weighted=False):
    """ Compute rank similarity of peak arrays of `pk1`, `pk2`
      where two jaccard simlilar peaks with approximate intesity ranks are considered similar
    @parms 
      d: number of digits to round up [NOT USED]
      delta: threshold for idential peaks [NEW]
      good_with: how near two ranks should be so that they can be considered as coincident
      weighted: whether to weigh higher/lower ranks
    """
    # npk1, npk2 = _peak(pk1, g=d), _peak(pk2, g=d)
    npk1, npk2 = pke.align(pk1, pk2, delta)
    common = npk1.keys() & npk2.keys()
    union  = npk1.keys() | npk2.keys()
    
    # n_common = # of co-occurrent "good" peaks  
    rpk1, rpk2 = pke.rank(npk1), pke.rank(npk2)
    good_peak = lambda x: abs(rpk1[x]-rpk2[x]) < good_with
    if weighted:
        weight = lambda r: 1/float(r) #  weighted by `reciprocal rank`: [1/n,1]
        good_common = [0.1+weight(rpk1[k])+weight(rpk2[k]) for k in common if good_peak(k)] 
        return sum(good_common)/len(union)
    else:
        good_common = [good_peak(k) for k in common] 
        return sum(good_common)/len(union)    

"""cf. https://en.wikipedia.org/wiki/Rank_correlation
"""