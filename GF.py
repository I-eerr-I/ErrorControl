class DPoly:
    def __init__(self, poly, mod = None):
        if type(poly) not in [list, DPoly]:
            if type(poly) != int:
                raise Exception("The first parameter must be type of list or int")
            else:
                poly = [0 for _ in range(poly)] + [1]
        elif type(poly) == DPoly:
            poly = poly.poly.copy()
        
        if mod:
            if type(mod) in [DPoly, list]:
                if type(mod) == list:
                    mod = DPoly(mod.copy())
                self.poly = (DPoly(poly.copy())/mod)[0]
            else:
                self.poly = (DPoly(poly.copy())/DPoly(mod))[0]
        else:
            self.poly = poly.copy()
        self.length = len(poly)
        self.deg    = self.get_deg()
        self.weight = self.get_weight(poly)
        self.update()

    def __rmul__(self, other):
        if type(other) in [int, float, bool]:
            return self.__mul__(other)
        else:
            other.__mul__(self)

    def __radd__(self, other):
        if type(other) in [int, float, bool]:
            return self.__add__(other)
        else:
            other.__add__(self)

    def copy(self):
        return DPoly(self.poly.copy())

    def get_deg(self):
        for i in range(len(self.poly)-1, 0, -1):
            if self[i] != 0:
                return i
        return 0    

    def get_weight(self, poly):
        w = 0
        for p in poly:
            if p != 0:
                w += 1
        return w

    def update(self):
        self.weight = self.get_weight(self.poly)
        self.length = len(self.poly)
        self.deg    = self.get_deg()
        for i in range(self.length):
            if type(self[i]) == GFEps and self[i].power == 0:
                self[i] = 1

    def __call__(self, value):
        to_add = []
        for i in range(self.length):
            if i == 0 and self[i] != 0:
                to_add.append(self[i])
            if self[i] != 0 and i != 0:
                to_add.append((value**i) * self[i])
        if self.deg == 0:
            return self[0]
        if len(to_add) == 1:
            return to_add[0]
        result = to_add[0] + to_add[1]
        for i in range(2,len(to_add)):
            result = result + to_add[i]            
        return result

    def __sub__(self, poly2):
        if type(poly2) == list:
            poly2 = DPoly(poly2.copy())
        return self + poly2

    def __rsub__(self, poly2):
        return poly2 + self

    def __add__(self, poly2, prnt = False):
        if type(poly2) not in [DPoly, list, tuple, set, dict]:
            result = self.poly.copy()
            result[0] = poly2 + result[0]
            return DPoly(result)

        if type(poly2) == list:
            poly2 = DPoly(poly2.copy())
        
        result = DPoly([0 for _ in range(max(len(self), len(poly2)))])
        for i in range(min(len(self), len(poly2))):
            if type(poly2[i]) == GFEps:
                result[i] = poly2[i].__add__(self[i])
            elif type(self[i]) == GFEps:
                result[i] = self[i].__add__(poly2[i])
            else:
                result[i] = self.poly[i]^poly2[i]
        if len(self) != len(poly2):
            if len(self) > len(poly2):
                for i in range(len(poly2), len(self)):
                    result[i] = self[i]
            else:
                for i in range(len(self), len(poly2)):
                    result[i] = poly2[i]
        result.update()
        return result

    def __mul__(self, poly2):
        if type(poly2) != DPoly and (poly2 != 0 or poly2 != 1):
            result = self.poly.copy()
            for i in range(len(self)):
                if type(poly2) != GFEps:
                    if self[i] != 1 and self[i] != 0:
                        result[i] = (poly2%2) * result[i]
                    elif self[i] == 1:
                        result[i] = poly2%2
                else:
                    if self[i] != 1 and self[i] != 0:
                        result[i] = (poly2) * result[i]
                    elif self[i] == 1:
                        result[i] = poly2
            return DPoly(result)
        if type(poly2) == list:
            poly2 = DPoly(poly2.copy())
        result  = DPoly([0 for _ in range(self.length+poly2.length-1)])
        if len(self) < len(poly2):
            return poly2.__mul__(self)
        
        if not self.weight or not poly2.weight:
            return result
        elif self.weight == 1 and self.poly[0] == 1:
            return DPoly(poly2.poly.copy())
        elif poly2.weight == 1 and poly2[0] == 1:
            return DPoly(self.poly.copy())

        sub_results = []
        for i in range(self.length):
            if i == 0 and self.poly[i] == 1:
                sub_results.append(poly2)
                continue
            elif self.poly[i] == 0:
                continue
            else:
                sub_result = DPoly([0 for _ in range(self.length + poly2.length-1)])
                for j in range(poly2.length):
                    if self.poly[i] * poly2[j] != 0:
                        sub_result[j+i] = self[i] * poly2[j]
                sub_results.append(sub_result)

        for i in range(len(sub_results)):
            result = sub_results[i] + result

        result.update()
        return result

    def __getitem__(self, index):
        return self.poly[index]

    def __setitem__(self, index, value):
        self.poly[index] = value

    def __repr__(self):
        return str(self.poly) + " DPoly"

    def equal_null(self):
        for pi in self.poly:
            if pi != 0:
                return False
        return True

    def __truediv__(self, poly2):
        if type(poly2) == list:
            poly2 = DPoly(poly2.copy())
        if type(poly2) != DPoly:
            return self / DPoly([poly2])
        if self.deg < poly2.deg:
            return [None, None]
        if poly2 == 1 or (poly2.weight == 1 and poly2[0] == 1):
            return [DPoly(self.poly.copy()), DPoly([0])]
        elif poly2 == 0 or (poly2.weight == 0):
            return [None, None]
        elif self.weight == 0:
            return [DPoly([0]),DPoly([0])]

        result   = DPoly([0])
        reminder = DPoly(self.poly.copy())
        while reminder.get_deg() >= poly2.get_deg():
            if type(reminder[reminder.deg]) != GFEps:
                if type(poly2[poly2.deg]) != GFEps:
                    poly = DPoly([0 for i in range(reminder.get_deg()-poly2.get_deg())] + [1])
                else:
                    poly = DPoly([0 for i in range(reminder.get_deg()-poly2.get_deg())] + [poly2[poly2.deg].to(poly2[poly2.deg].field(0))])
            else:
                if type(poly2[poly2.deg]) != GFEps:
                    poly = DPoly([0 for i in range(reminder.get_deg()-poly2.get_deg())] + [reminder[reminder.deg]])
                else:
                    poly = DPoly([0 for i in range(reminder.get_deg()-poly2.get_deg())] + [poly2[poly2.deg].to(reminder[reminder.deg])])
            result = poly + result
            poly = poly * poly2
            reminder = reminder + poly
            if reminder.equal_null():
                break
        result.update()
        reminder.update()
        return [result, reminder]

    def mul_mod(self, poly2, mod):
        return DPoly((self * poly2).poly, mod)

    def __str__(self):
        if self.weight == 0:
            return "0"
        if self.weight == 1 and self[0] != 0:
            return str(self[0])
        poly_str = str(self[0]) + " + " if self[0] != 0 else ""
        counter = 1
        for i in range(1, len(self)):
            if self[i] == 0:
                continue
            if self[i] != 1:
                poly_str += str(self[i])+" * x^"+str(i)
            else:
                poly_str += "x^"+str(i)
            if self.weight != counter:
                poly_str += " + "
            counter += 1
        if poly_str[-2] == "+":
            poly_str = poly_str[:-2]
        return poly_str

    def __len__(self):
        return self.length

    def der(self, times = 1):
        poly_ = self[1]
        for j in range(2, self.deg+1):
            poly_ += (j%2)*self[j]*DPoly(j-1)
        if times - 1 > 0:
            return self.der(times-1)
        return poly_

class GFEps:

    def __init__(self, poly, field, power):
        self.field  = field
        self.eps    = poly.copy()
        self.deg    = field.deg
        self.power  = power 

    def to(self, eps):
        return self.field(eps.power - self.power) 

    def __getitem__(self, index):
        return self.eps[index]

    def __rmul__(self, other):
        if type(other) in [int, float, bool]:
            return self.__mul__(other)
        else:
            other.__mul__(self)

    def __radd__(self, other):
        if type(other) in [int, float, bool]:
            return self.__add__(other)
        else:
            other.__add__(self)

    def __abs__(self):
        return self

    def __len__(self):
        return self.deg

    def __add__(self, other):
        result = [0 for _ in range(self.deg)]
        if type(other) == GFEps:
            if self.power == other.power:
                return 0
            for i in range(self.deg):
                result[i] = self[i]^other[i]
            return GFEps(result, self.field, self.field.power_of(result))
        elif other == 0:
            return GFEps(self.eps.copy(), self.field, self.power)
        elif other == 1:
            return self + self.field(0)
        elif type(other) == DPoly:
            return other.__add__(self)
        else:
            try:
                power = self.field.power_of(other, base=10)
                return self + self.field(power) 
            except ValueError:
                return None

    def __mul__(self, other):
        if type(other) == DPoly:
            return other.__mul__(self)
        elif type(other) != GFEps:
            if other == 1:
                return GFEps(self.eps.copy(), self.field, self.power)
            elif other == 0:
                return 0
            elif other <= -1:
                return self * abs(other)
            else:
                try:
                    power = self.field.power_of(other, base=10)
                    return self * self.field(power) 
                except ValueError:
                    return None
        power = (self.power + other.power) % (2**self.deg - 1) 
        return self.field(power)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "eps^"+str(self.power)

    def __pow__(self, power):
        power *= self.power
        power %= 2**self.deg - 1
        return self.field(power)

    def __eq__(self, other):
        if type(other) not in [list, tuple, set, dict, GFEps, DPoly]:
            try:
                return self == self.field(self.field.power_of(other))
            except:
                return False
        for i in range(len(self)):
            if self[i] != other[i]:
                return False
        return True

polynom = {1:  DPoly([1,1]),
           2:  DPoly([1,1,1]),
           3:  DPoly([1,1,0,1]),
           4:  DPoly([1,1,0,0,1]),
           5:  DPoly([1,0,1,0,0,1]),
           6:  DPoly([1,1,0,0,0,0,1]),
           7:  DPoly([1,1,0,0,0,0,0,1]),
           8:  DPoly([1,0,1,1,1,0,0,0,1]),
           9:  DPoly([1,0,0,0,1,0,0,0,0,1]),
           10: DPoly([1,0,0,1,0,0,0,0,0,0,1]),
           11: DPoly([1,0,1,0,0,0,0,0,0,0,0,1]),
           12: DPoly([1,1,0,0,1,0,1,0,0,0,0,0,1]),
           13: DPoly([1,1,0,1,1,0,0,0,0,0,0,0,0,1]),
           14: DPoly([1,1,0,1,0,1,0,0,0,0,0,0,0,0,1]),
           15: DPoly([1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1]),
           16: DPoly([1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1]),
           17: DPoly([1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1]),
           18: DPoly([1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1]),
           19: DPoly([1,1,1,0,0,1])+DPoly(19),
           20: DPoly([1,0,0,1])+DPoly(20),
           21: DPoly([1,0,1])+DPoly(21),
           22: DPoly([1,1])+DPoly(22)}

class GF:
    
    def __init__(self, poly, deg = None):
        if deg:
            self.deg     = deg
        else:
            if type(poly) != DPoly:
                if type(poly) == int:
                    poly = polynom[poly]
                else:
                    poly = DPoly(poly.copy())
            self.deg = poly.deg
        self.poly    = poly.copy()
        self.length  = 2**self.deg - 1
        self.field10 = []
        self.field   = []
        self.build_field()
    
    def __call__(self, index, base = 2):
        if base == 10:  return self.field10[index]
        elif base == 2:
            if index > len(self.field)-1:
                index = index % (2**self.deg - 1)
            return self.field[index]
    
    def build_field(self):
        zeros = [0 for _ in range(self.deg-1)]
        eps1  = zeros.copy()
        eps1.insert(0,1)
        eps2  = zeros.copy()
        eps2.insert(1,1)
        eps3  = zeros.copy()
        eps3.insert(2,1)
        self.field.append(GFEps(eps1, self, 0))
        self.field.append(GFEps(eps2, self, 1))
        self.field.append(GFEps(eps3, self, 2))
        
        for i in range(2, 2**self.deg - 2):
            self.field.append(self.mul_poly(self.field[1], self.field[i], i+1))
        for i in range(len(self.field)):
            str_res = ""
            for x in self.field[i]:
                str_res += str(x)
            self.field10.append(int(str_res, 2))

    def pow_ind(self, index, power):
        return (index*power) % (self.deg**2 - 1)

    def mul_poly(self, poly1, poly2, power = None):
        result  = [0 for _ in range(self.deg)]
        weight1 = self.get_weight(poly1)
        weight2 = self.get_weight(poly2)

        if not weight1 or not weight2:
            return GFEps(result, self, power)
        elif weight1 == 1 and poly1[0] == 1:
            return GFEps(poly2, self, power)
        elif weight2 == 1 and poly2[0] == 1:
            return GFEps(poly1, self, power)

        sub_results = []
        for i in range(self.deg):
            if i == 0 and poly1[i] == 1:
                sub_results.append(poly2)
                continue
            elif poly1[i] == 0:
                continue
            else:
                mod_result = []
                sub_result = [0 for _ in range(self.deg)]
                for j in range(self.deg):
                    if poly1[i] * poly2[j] == 1:
                        if j == self.deg-1:
                            mod_result = self.poly[0:self.deg]
                        else:
                            sub_result[j+i] = 1
                if mod_result:
                    sub_results.append(self.add_poly(mod_result, sub_result, i))
                else:
                    sub_results.append(sub_result)
        for sub_res in sub_results:
            result = self.add_poly(sub_res, result, i).eps
        return GFEps(result, self, power)
        

    def add_poly(self, poly1, poly2, power = None):
        result = [0 for _ in range(self.deg)]
        for i in range(self.deg):
            result[i] = poly1[i]^poly2[i]
        return GFEps(result, self, power)
                
    def power_of(self, value, base = 2):
        if value == [0 for _ in range(self.deg)]: return -1
        if base == 2:
            if type(value) == GFEps:
                for i in range(len(self.field)):
                    if self.field[i].eps == value.eps:
                        return i
            else:
                for i in range(len(self.field)):
                    if self.field[i].eps == value:
                        return i
        if base == 10: return self.field10.index(value)

    def get_weight(self, poly):
        w = 0
        for p in poly:
            if p == 1:
                w += 1
        return w


if __name__=='__main__':
    a = DPoly([0,1,1])
    b = DPoly([0,0,1])
    print(a*b)
    gf = GF(DPoly([1,1,0,0,1]))
    for eps in gf.field: print(eps.eps)
    c = a * gf(1)
    d = c + gf(2)
    print(d*b)
    a = DPoly([gf(1), 1])
    r1, r2 = b/a
    
