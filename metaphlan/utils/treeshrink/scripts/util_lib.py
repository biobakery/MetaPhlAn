from math import log, sqrt

def minVar_outliers(L,max_out_num,min_gvf=0.2):
# perform minVar bisect multiple times as long as 
# the lower half (i.e. L1) can still be bisect with 
# gvf at least the specified min_gvf; if the first 
# bisection has gvf < min_gvf, then return the max value
    L.sort()
    n = len(L)
    prev_cutoff, prev_cutoff_idx, prev_gvf = None,None,None
    cutoff, cutoff_idx, gvf = minVar_bisect(L)
    while gvf is not None and gvf >= min_gvf and n-cutoff_idx-1 <= max_out_num:
        prev_cutoff, prev_cutoff_idx, prev_gvf = cutoff, cutoff_idx, gvf
        cutoff, cutoff_idx, gvf = minVar_bisect(L[:cutoff_idx+1])

    if prev_gvf is None:
        cutoff_idx = len(L)-1
        cutoff = L[cutoff_idx]
        gvf = 0
    else:
        cutoff, cutoff_idx, gvf = prev_cutoff, prev_cutoff_idx, prev_gvf    

    return cutoff, cutoff_idx, gvf    

def minVar_bisect(L):
# cut a list L into two sets L1, L2 such that
# the var(L1) + var(L2) is minimized (Jenks natural 2-break)
# return (max(L1)+min(L2))/2 as the cut-off
# NOTE: assume that list L is already sorted

    sum_left = 0
    sumsq_left = 0
    sum_right = 0
    sumsq_right = 0

    for x in L:
        sum_left += x
        sumsq_left += x*x

    n = float(len(L))

    var0 = sumsq_left/n - (sum_left/n)**2
    minVar = var0
    cutoff = L[-1]
    
    k = len(L)-1
    cutoff_idx = k

    while k > 0:
        x = L[k]    
        sum_left -= x
        sumsq_left -= x*x
        sum_right += x
        sumsq_right += x*x
     
        var_left = sumsq_left/float(k) - (sum_left/float(k))**2
        var_right = sumsq_right/(n-k) - (sum_right/(n-k))**2
        var = var_left+var_right
        if var < minVar:
            minVar = var
            cutoff = (L[int(k-1)]+L[int(k)])/2
            cutoff_idx = k-1
        k -= 1    
    gvf = (var0-minVar)/var0 if var0 else None 
    return cutoff, cutoff_idx, gvf

def minCV_bisect(L):
# cut a list L into two sets L1, L2 such that
# the sum of the coefficient of variation of 
# L1 and L2 is minimized 
# return (max(L1)+min(L2))/2 as the cut-off
# NOTE: assume that list L is already sorted

    sum_left = 0
    sumsq_left = 0
    sum_right = 0
    sumsq_right = 0

    for x in L:
        sum_left += x
        sumsq_left += x*x

    n = float(len(L))

    var0 = sumsq_left/n - (sum_left/n)**2
    cv0 = sqrt(var0)/(sum_left/n)
    minCV = cv0
    cutoff = L[-1]
    
    k = len(L)-1
    cutoff_idx = k

    while k > 0:
        x = L[k]    
        sum_left -= x
        sumsq_left -= x*x
        sum_right += x
        sumsq_right += x*x
     
        var_left = sumsq_left/float(k) - (sum_left/float(k))**2
        cv_left = sqrt(abs(var_left))/(sum_left/float(k))
        var_right = sumsq_right/(n-k) - (sum_right/(n-k))**2
        cv_right = sqrt(abs(var_right))/(sum_right/(n-k))
        cv = cv_left + cv_right
        if cv < minCV:
            minCV = cv
            cutoff = (L[int(k-1)]+L[int(k)])/2
            cutoff_idx = k-1
        k -= 1    
    #gvf = (var0-minVar)/var0 if var0 else None 
    return cutoff, cutoff_idx
