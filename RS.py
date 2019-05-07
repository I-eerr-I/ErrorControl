from GF import GF, DPoly
from NIAlgorithm import find_syndromes, BMA, PGCA, EuA, ChenA, ForniA
from math import sqrt
from random import random

class RS:
    def __init__(self, n, t = 2, b = 1):
        self.n = n
        self.t = t
        self.b = b
        self.k = n - 2*t
        self.field = GF(int(sqrt(n+1)))
        self.g = DPoly([self.field(b), 1])
        for j in range(b+1, b+2*t):
            self.g *= DPoly([self.field(j), 1])
        self.h = ((DPoly([1]) + DPoly(15))/self.g)[0]

    def encode(self, u, sys = True):
        if sys:
            shifted_u = u*DPoly(self.n - self.k)
            res, r = shifted_u/self.g
            return shifted_u + r
        else:
            return self.g*u

    def decode(self, r, sys = True, debug = False):
        if sys:
            S = find_syndromes(r, self.field, 2*self.t)
            Lambda_BMA  = BMA(S, self.t, self.field(0))
            Lambda_PGCA = PGCA(S, self.t)
            Lambda_EuA  = EuA(S, self.t)
            count = 0
            bma_pgca = False
            bma_eua  = False
            pgca_eua = False
            if debug:
                print("BMA", Lambda_BMA)
                print("PGCA", Lambda_PGCA)
                print("EuA", Lambda_EuA)
            try:
                if type(Lambda_PGCA) == DPoly and type(Lambda_BMA) == DPoly and type(Lambda_EuA) == DPoly:
                    try:
                        assert Lambda_BMA.poly  == Lambda_PGCA.poly
                        bma_pgca = True
                    except:
                        pass
                    try:
                        assert Lambda_PGCA.poly == Lambda_EuA.poly
                        pgca_eua = True
                    except:
                        pass
                    try:
                        assert Lambda_EuA.poly  == Lambda_BMA.poly
                        bma_eua = True
                    except:
                        pass
                    
                    if Lambda_BMA.poly == [1]:
                        return (((r/DPoly(self.n-self.k))[0]*self.g)/(self.g + DPoly([1])))[0]
                else:
                    try:
                        assert Lambda_BMA  == Lambda_PGCA
                        bma_pgca = True
                    except: pass
                    try:
                        assert Lambda_PGCA == Lambda_EuA
                        pgca_eua = True
                    except: pass
                    try:
                        assert Lambda_EuA  == Lambda_BMA
                        bma_eua  = True
                    except: pass
            except Exception as e:
                if debug:
                    print(e)
                return (((r/DPoly(self.n-self.k))[0]*self.g)/(self.g + DPoly([1])))[0]
            if bma_pgca or bma_eua:
                Lambda = Lambda_BMA
            elif pgca_eua:
                Lambda = Lambda_PGCA
            else:
                return (((r/DPoly(self.n-self.k))[0]*self.g)/(self.g + DPoly([1])))[0]
            S_x = DPoly([S[i] for i in S])
            Omega = ((S_x*Lambda)/DPoly(2*self.t))[1]
            roots = ChenA(Lambda, self.field)
            errors = ForniA(roots, Omega, Lambda, self.t, self.b)
            error = DPoly(self.n)
            error[self.n] = 0
            error.update()
            for i in range(len(roots)):
                error[roots[i].power] = errors[i]
                error.update()
            v = r + error
            v.update()
            return (((v/DPoly(self.n-self.k))[0]*self.g)/(self.g + DPoly([1])))[0]
        else:
            return None

if __name__=="__main__":
    rs = RS(15, 3)
    gf = GF(4)
    p  = 0.05
    print(rs.k)
    u = DPoly([gf(i) for i in range(rs.k)])
    print(u)
    v = rs.encode(u)
    print(v)
    count = 0
    for i in range(len(v)):
        if random() < p and count < 2:
            v[i] = v[i]**-1
        count += 1
    print(rs.decode(v, debug=True))
    
