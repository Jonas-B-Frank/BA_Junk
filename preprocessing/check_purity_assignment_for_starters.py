from scipy.stats import binom
import math

purity = 0.9999
x = 0
y = 2
tcn_tumor = 2
k = 5
n = 10
pl_1 = 'nan'
pl_2 = 'nan'

if purity == 1 and ((x == 0) or (y == 0)) and n > 0:
    if x == 0:
        ref = x
        alt = y
    else:
        alt = x
        ref = y
    quality_score = 0

else:
    # case 1: x,y (ref, alt)
    vaf_1 = ((y * purity) + (1 - purity)) / ((tcn_tumor * purity) + (2 * (1 - purity)))
    prob_1 = binom.pmf(k=k, n=n, p=vaf_1)
    
    # case 2: y,x (ref, alt)
    vaf_2 = ((x * purity) + (1 - purity)) / ((tcn_tumor * purity) + (2 * (1 - purity)))
    prob_2 = binom.pmf(k=k, n=n, p=vaf_2)

    # calculating raw phred-scaled likelihoods
    pl_1 = -10 * math.log10(prob_1)
    pl_2 = -10 * math.log10(prob_2)

    if prob_1 > prob_2:
        ref = x
        alt = y
    else:
        ref = y
        alt = x

print(ref, alt, pl_1, ' ', pl_2)