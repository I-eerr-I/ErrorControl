from GF import DPoly, GF
from NIAlgorithm import BMA, PGCA, ChenA, EuA, find_syndromes
from math import sqrt

class BCH:

    def __init__(self, n, t = 2):
        self.n = n
        self.field = GF(int(sqrt(n+1)))
        self.t = t
        self.get_params()

    def get_params(self):
        minpolies = self.get_minpolies()
        self.g = minpolies[1]
        for i in range(2, 2*self.t+1):
            if i in minpolies:
                self.g = self.g * minpolies[i]
        self.k = self.n - self.g.deg

    def get_minpolies(self):
        roots = self.get_minpolies_roots()
        minpolies = {}
        for power in roots:
            minpolies[power] = DPoly([0,1]) + DPoly([roots[power][0]])
            for i in range(1, len(roots[power])):
                minpolies[power] = minpolies[power] * (DPoly([0,1]) + DPoly([roots[power][i]]))
        return minpolies
                
    
    def get_minpolies_roots(self):
        roots = {}
        roots[0] = [self.field(0)]
        found = []
        found.append(0)
        for i in range(1, self.n):
            if i not in found:
                roots[i] = [self.field(i)]
                found.append(i)
                cycle = [i]
                power = 1
                q = 0
                while True:
                    q = (2**power)*i % self.n
                    if q not in cycle:
                        cycle.append(q)
                        found.append(q)
                        roots[i].append(self.field(q))
                        power += 1
                    else:
                        break
        return roots

    def encode(self, inf):
        return self.g * inf

    def decode(self, encoded, debug = False):
        S = find_syndromes(encoded, self.field, 2*self.t)
        BMA_r  = BMA(S, self.t, self.field(0), debug)
        PGCA_r = PGCA(S, self.t, debug)
        EuA_r  = EuA(S, self.t, debug)
        if type(BMA_r) == DPoly and type(PGCA_r) == DPoly and type(EuA_r) == DPoly: 
            if BMA_r.poly == [1] and PGCA_r.poly == [1] and EuA_r == [1]:
                return (encoded / self.g)[0]
        for i in range(min(len(BMA_r), len(PGCA_r), len(EuA_r))):
            if BMA_r[i] != PGCA_r[i] != EuA_r[i]:
                return None

        Lambda = EuA_r
        roots = ChenA(Lambda, self.field)
        error = DPoly([0 for _ in range(self.n)])
        for root in roots:
            error[root.power] = 1
        encoded = error + encoded
        return (encoded / self.g)[0]
            

if __name__ == "__main__":
    bch = BCH(15, 3)
    print(bch.k)
    u   = DPoly([0,0,1,0,1])
    r   = bch.encode(u)
    v   = r + DPoly([0,0,0,1,0,0,0,0,0,1])
    u_d = bch.decode(v)
    print(u)
    print(r)
    print(v)
    print(u_d)
    assert u.poly == u_d.poly
