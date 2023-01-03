# 類似度計算
import math


def similarity(pk1, pk2, method='jaccard'):
    npk1, npk2 = _peak(pk1), _peak(pk2)
    common =  npk1.keys() & npk2.keys()
    diff1 = npk1.keys() - npk2.keys()
    diff2 = npk2.keys() - npk1.keys()
    diff  = diff1 | diff2
    n_common, n_diff1, n_diff2 = len(common), len(diff1), len(diff2)
    n = n_common + n_diff1 + n_diff2
    
    """cf. https://en.wikipedia.org/wiki/Rank_correlation
    """
    # spearman's rank correlation coefficient
    if method == 'spearman':       
        r1, r2 = _rank(npk1), _rank(npk2)
        rnk = {k: (r1[k]-r2[k])**2 for k in common}
        return 1 - sum(rnk.values())/float(n_common*(n_common-1))
    
    # kendall's rank correlation coefficient
    if method == 'kendall':
        return n_common / float(n)
    
    # defalut jaccard similarity
    return n_common / float(n)

# m/zの値は正確ではないため、下記の関数で結果を丸める
def _fmt(f, d=1):
    """ 実数を整数に丸める。d:桁数 d桁以内四捨五入 
       例：fmt(3456.78,d=1) == 3460, fmt(3456.78,d=2) == 3500
    """
    return round(int(f), -d)

def _peak(pk, d=1):
    """ Peaksのm/z値を丸め、辞書にして返す
    """
    npk =  ([_fmt(x,d=d) for x in pk[0]],pk[1])
    return dict(zip(npk[0], npk[1]))

def _rank(pk):
    """returns a new dict with rank for each m/z value 
    """
    return {key: rank for rank, key in enumerate(sorted(pk, key=pk.get, reverse=True), 1)}

