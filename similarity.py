# 類似度計算
import math

def jaccard_similarity(pk1, pk2):
    npk1, npk2 = _peak(pk1, d=1), _peak(pk2, d=1)
    common = npk1.keys() & npk2.keys()
    union  = npk1.keys() | npk2.keys()

    # n_common: # of co-occurent peaks
    n_common, n_union = len(common), len(union)

    return n_common / float(n_union)

def rank_similarity(pk1, pk2, good_with=5):
    npk1, npk2 = _peak(pk1, d=1), _peak(pk2, d=1)
    common = npk1.keys() & npk2.keys()
    union  = npk1.keys() | npk2.keys()
    
    # n_common = # of co-occurrent "good peak"s whose ranks in 2 series are close enough 
    rpk1, rpk2 = _rank(npk1), _rank(npk2)
    good_peak = lambda x: abs(rpk1[x]-rpk2[x]) < good_with
    good_common = [good_peak(k) for k in common]
    n_common, n_union = sum(good_common), len(union)
    
    return n_common / float(n_union)
    
"""cf. https://en.wikipedia.org/wiki/Rank_correlation
"""

def _peak(pk, d=1):
    """ Peaksのm/z値を丸め、辞書にして返す
       実数を整数に丸める。d:桁数 d桁以内四捨五入 
       例：_fmt(3456.78,d=1) == 3460, _fmt(3456.78,d=2) == 3500
    """
    npk =  ([round(int(x), -d) for x in pk[0]],pk[1])
    return dict(zip(npk[0], npk[1]))

def _rank(pk):
    """returns a new dict with rank for each m/z value 
    """
    return {key: rank for rank, key in enumerate(sorted(pk, key=pk.get, reverse=True), 1)}

