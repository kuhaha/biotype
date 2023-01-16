# 類似度計算
import math

def jaccard_similarity(pk1, pk2, d=1):
    npk1, npk2 = _peak(pk1, g=d), _peak(pk2, g=d)
    common = npk1.keys() & npk2.keys()
    union  = npk1.keys() | npk2.keys()

    # n_common: # of co-occurent peaks
    n_common, n_union = len(common), len(union)
    return n_common / float(n_union)

def rank_similarity(pk1, pk2, d=1, good_with=2, weighted=False):
    npk1, npk2 = _peak(pk1, g=d), _peak(pk2, g=d)
    common = npk1.keys() & npk2.keys()
    union  = npk1.keys() | npk2.keys()
    
    # n_common = # of co-occurrent "good" peaks  
    rpk1, rpk2 = _rank(npk1), _rank(npk2)
    good_peak = lambda x: abs(rpk1[x]-rpk2[x]) < good_with
    if weighted:
        weight = lambda r: 1/float(r) #  weighted by `reciprocal rank`: [1/n,1]
        good_common = [0.1+weight(rpk1[k])+weight(rpk2[k]) for k in common if good_peak(k)] 
#         n1,n2 = len(rpk1), len(rpk2)
#         weight = lambda r, n: 1-(r-1)/float(n) # weighted by `rank sum`: [1/n,1]
#         good_common = [0.5+weight(rpk1[k],n1)+weight(rpk2[k],n2) for k in common if good_peak(k)] 
        return sum(good_common)/len(union)
    else:
        good_common = [good_peak(k) for k in common] 
        return sum(good_common)/len(union)    

"""cf. https://en.wikipedia.org/wiki/Rank_correlation
"""

def _peak(pk, g=1):
    """ Peaksのm/z値を丸め、辞書にして返す
       実数を整数に丸める。g:桁数 g桁以内四捨五入 
       例：_fmt(3456.78, g=1) == 3460, _fmt(3456.78, g=2) == 3500
    """
    npk =  ([round(int(x), -g) for x in pk[0]],pk[1])
    return dict(zip(npk[0], npk[1]))

def _rank(pk):
    """returns a new dict with rank for each m/z value 
    """
    return {key: rank for rank, key in enumerate(sorted(pk, key=pk.get, reverse=True), 1)}

