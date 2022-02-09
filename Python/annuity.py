from math import floor
from numpy.matlib import repmat


def schedule(coupon, principal, years, paymentsPerYear):
    n = paymentsPerYear
    T = floor(years*n)/n
    #residual = years - T
    R_tilde = coupon / n
    Y = R_tilde/(1-(1+R_tilde)**(-n*T))*principal
    F_t = (1 + R_tilde). ^ i * principal - Y * ((1 + R_tilde). ^ i - 1). / R_tilde

