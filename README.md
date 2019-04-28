# GF

This script includes three classes: **GF, GFEps, DPoly.**
---
### GF class 
*is realization of extended Galois Field on binary base - GF(2^m). Using GF(degree = m, polynomial = p(x)) you create new Galois field on binary base with the power "m" and generating polynomal - p(x), where p(x) is the list with 1 and 0 from lower to higher power. Your created G field will contain all elements of type GFEps, that you can multiply, divide, add and etc with each other. For example:*

**Create G field**
```
1. p_x = [1,1,0,0,1] # p(x) = 1+x+x^4
2. m   = 4
3. gf  = GF(m, p_x)
```

Now you can access all elements by calling gf with the power of elements you want to get:

**Access elements**
```
4. gf(0)  # eps^0 = 1 - the element of GF(2^4) with 0 power
5. gf(1)  # eps^1
6. gf(15) # eps^15 = esp^0 = 1 mod(2^4 - 1)
7. gf(24) # eps^24 = eps^9 mod(2^4 - 1)
8. gf(-1) # eps^14
```
and etc

GF class has 4 methods to build field, that are not very usefull for simple work with created field. But there are 3 methods that you can use by your own:

**pow_ind(index, power) - get the power of raised element's degree with index = index by power = power**
```
9.  gf.pow_ind(3, 4)  # 3*4 = 12, 12 % 2^4 - 1 = 12
10. gf.pow_ind(10, 2) # 10*2 = 20, 20 % 2^4 - 1 = 5
```

**power_of(value, base = 2) - get the power of element = value using the base of 2 or 10**
```
11. gf.power_of([1,0,0,0]) # return is 0
12. gf.power_of([0,1,0,0]) # return is 1
```

**subs_value(poly, val_power) - substitude the element of power val_power in the polinomial poly** 
*the poly should be with type of DPoly*
```
13. gf.subs_value(DPoly([0,1]), 1) # DPoly([1]) = x = f(x), subs_value(f(x), 1) => f(gf(1)) = f(eps^1) = eps^1
14. gf.subs_value(DPoly([0,1,1]), 2) # DPoly([0,1,1]) = x+x^2 = h(x), subs_value(h(x), 2) => h(gf(2)) = h(eps^2) = eps^2 + eps^4
```

