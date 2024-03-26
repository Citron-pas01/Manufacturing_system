import numpy as np
import scipy.stats as st
import math

coefficients = [ 4.41738119e-09, 1.79200966e-07, 3.01634229e-06,
                2.63537452e-05, 1.12381749e-04,   5.71289020e-06,
               -2.64198510e-03,  -1.59986142e-02, -5.60399292e-02,
                -1.48968884e-01, -3.68776346e-01, -1.22551895e+00,
               -8.99375602e-01]
def inverse_standard_loss(L):
    x = math.log(L)
    z = np.polyval(coefficients, x)
    return z


'''
C = 1.5
I = 0.28
A = 100
P = 12.8
d = 280 # weeek
sigma_0 = np.sqrt(77*77)
l = 5 # week
'''

C = 18.8
I = 0.4
A = 75
P = 400
d = 38 # weeek
sigma_0 = np.sqrt(130)
l = 3 # week
beta = 0.99

D = d*52
h = I*C
mu = d*l
sigma = np.sqrt(l)*sigma_0


Q_EOQ = np.sqrt(2*A*D/h)
Q = math.ceil(Q_EOQ)
print('Q',Q)
Q_n_1 = Q-3

while Q-Q_n_1 >2:
    print('again')
    F_r = 1-Q*h/(P*D)
    Q_n_1 = Q
    Z = st.norm.ppf(F_r)
    #print(F_r,Z)
    R = math.ceil(mu + sigma*Z)
    L_z = st.norm.pdf(Z) - Z*(1-st.norm.cdf(Z))
    n_r = sigma*L_z
    #print('n_r', n_r)
    Q_n = np.sqrt(2*D*(A+P*n_r)/h)
    Q = math.ceil(Q_n)
    #print('(Q,R) is ({0},{1})'.format(Q,R))

print('Final (Q,R) is ({0},{1})'.format(Q,R))

if beta > 0:
    #Initialize Intermediate variables 
    Qinter = 0
    Rinter = 0
    Rafter = 1
    #Estimate Q with EOQ
    Qafter = Q_EOQ
    while (round(Qinter)!= round(Qafter)
    and round(Rinter)!= round(Rafter)): #Set n(R)
        n_R = (1-beta)*Qafter #Calculate L(z)
        L_z = n_R/sigma
        #Find Z-score
        Zn = inverse_standard_loss(L_z) #Find 1-F(z)
        One_F = 1-st.norm.cdf(Zn) #Store previous R and Q
        Rinter = Rafter
        Qinter = Qafter
        #Find new R
        Rafter = sigma * Zn + mu #Find new Q
        Qafter = n_R/One_F + np.sqrt(2*A*D/h+(n_R/One_F)*(n_R/One_F))
#Display Results
    print ("R optimal = " + str(round(Rafter)))
    print ("Q optimal = " + str(round(Qafter)))