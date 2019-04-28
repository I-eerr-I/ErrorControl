# GF
Python script to work with the Galois field

This script includes three classes: GFEps, GF, DPoly.
GF class is realization of extended Galois Field on binary base - GF(2^m). Using GF(degree = m, polynomial = p(x)) you create new Galois field on binary base with the power "m" and generating polynomal - p(x), where p(x) is the list with 1 and 0 from lower to higher power. Your created G field will contain all elements of type GFEps, that you can multiply, divide, add and etc with each other. For example:
p_x = [1,1,0,0,1] # p(x) = 1+x+x^4
m   = 4
gf  = GF(m, p_x)

Now you can access all elements by calling gf with the power of elements you want to get:

gf(0)  # eps^0 = 1 - the element of GF(2^4) with 0 power
gf(1)  # eps^1
gf(15) # eps^15 = esp^0 = 1 mod(2^4 - 1)
gf(24) # eps^24 = eps^9 mod(2^4 - 1)
gf(-1) # eps^14
and etc
