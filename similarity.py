# Compute similarity between peak arrays, where a peak is in the form <m/z, intensity>

import math

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
    npk1, npk2 = align(pk1, pk2, delta)
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
    npk1, npk2 = align(pk1, pk2, delta)
    common = npk1.keys() & npk2.keys()
    union  = npk1.keys() | npk2.keys()
    
    # n_common = # of co-occurrent "good" peaks  
    rpk1, rpk2 = _rank(npk1), _rank(npk2)
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

def _peak(pk, g=1): 
    """ Peaksのm/z値を丸め、辞書にして返す  !!! NOT USED, changed to align() 
      実数を整数に丸める。g:桁数 g桁以内四捨五入 
      例：_peak(3456.78, g=1) == 3460, _peak(3456.78, g=2) == 3500
    """
    npk =  ([round(int(x), -g) for x in pk[0]],pk[1])
    return dict(zip(npk[0], npk[1]))

def _rank(pk):
    """returns a new dict with rank for each m/z value 
    """
    return {key: rank for rank, key in enumerate(sorted(pk, key=pk.get, reverse=True), 1)}

def align(pk1, pk2, delta = 3):
    """ align two peak lists, assign a pair of m/z values in different lists
      an identical id if they are close enough (within `delta`)
      @params: 
      @returns: dictionary for each peak list, dict item {m/z: peak_id}
    """
    n1, n2 = len(pk1[0]), len(pk2[0])
    _pk1 = list(zip(pk1[0], [0]*n1)) # a list of tuples (m/z, 0) 
    _pk2 = list(zip(pk2[0], [1]*n2)) # a list of tuples (m/z, 1)
    _pk1.extend(_pk2)
    npk_sorted = sorted(_pk1, key=lambda x: x[0])
    n = len(npk_sorted)

    peak_id = 1
    npk = [dict(),  dict()]
    for i in range(n):
        m1, k1 = npk_sorted[i]
        if i < n - 1:
            m2, k2 = npk_sorted[i+1]
            if k1 != k2 and abs(m1-m2) <= delta:
                npk[k1][m1] = peak_id
                npk[k2][m2] = peak_id
                peak_id = peak_id + 1
        #       print(f"{k1}:{m1}, {k2}:{m2},delta={abs(m1-m2):.2f}")
        if m1 not in npk[k1].keys():
            npk[k1][m1] = peak_id
            peak_id = peak_id + 1
            
    #   print(f"{m1}: {npk[k1][m1]}")
    npk0_swap = {v: k for k, v in npk[0].items()}
    npk1_swap = {v: k for k, v in npk[1].items()}

    return npk0_swap, npk1_swap
