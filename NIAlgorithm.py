from DMatrix import DMatrix
from GF import DPoly, GF, GFEps

#Piterson-Gorenshtein-Zirler Algorithm for key equation
def PGCA(S, t, debug = False, return_data = False):
    Q   = 0
    det = 0
    Ms  = []
    for q in range(t, 0, -1):
        Ms = DMatrix([[S[i] for i in range(j, q+j)] for j in range(1, q+1)])
        if debug:
            print("Возьмем q = "+str(q)+" и составим матрицу синдромов Ms:")
            print(Ms)
        det = Ms.det()
        if debug:
            print("Найдем определить матрицы Ms:")
            print("det_"+str(q)+" = "+str(det))
        if det == 0:
            if debug:
                print("Так как det = 0, продолжим алгоритм, уменьшив q на единицу")
            continue
        else:
            if debug:
                print("Так как det != 0, то количество ошибок равно q = "+str(q))
                print("Переходим к решению ключевого уравнения")
            Q = q
            break
    if Q == 0:
        if debug:
            print("В принятой комбинации больше t ошибок или ошибок нет вообще!")
        return DPoly([1])
    Mi = Ms.inverse()
    M2 = DMatrix([[S[i] for i in range(q+1, 2*q+1)]])
    Mr = M2*Mi
    Lambda = Mr.matrix[0]
    Lambda.reverse()
    Lambda = DPoly([1] + [Lambda[i] for i in range(len(Lambda))])
    if debug:
        print("Матрица синдромов для q = "+str(q)+" равна Ms = ")
        print(Ms)
        print("Обратная ей матрица Mi = Ms^-1 = ")
        print(Mi)
        print("Матрица многочлена значений ошибок M2 = ")
        print(M2)
        print("Ключевое уравнение в матричном виде равно M2 = Ms*Ml")
        print("где Ml - матрица полинома локаторов ошибок. Таким образом")
        print("Ml = M2*Mi = ")
        print(Mr)
        print("От сюда получаем полином локаторов ошибок:")
        print("Lambda(x) = "+str(Lambda))
    data = {"Ms": Ms, "Mi": Mi, "M2": M2, "Lambda": Lambda}
    if return_data:
        return (data, Lambda)
    return Lambda

#Euclid algorithm for key equation
def EuA(S, t, debug = False, return_data = False):
    if type(S) != DPoly:
        if type(S) == list:
            S = DPoly(S)
        else:
            S = DPoly([1] + [S[i] for i in range(1, 2*t+1)])
    q = {0: DPoly([0 for i in range(2*t+1)] + [1]), 1: S}
    b = {0: DPoly([0]),  1: DPoly([1])}
    if debug:
        print("Начальные значения:")
        print("q[0] = x^(2t+1) = "+str(q[0]))
        print("q[1] = S(x) = "+str(q[1]))
        print("b[0] = 0")
        print("b[1] = 1")
    j = 1
    f = {}
    data = {"q": q, "b": b, "f": f, "j": j}
    while q[j].deg > t:
        j += 1
        data["j"] = j
        [f[j], q[j]] = q[j-2]/q[j-1]
        assert (0 <= q[j].deg and q[j].deg < q[j-1].deg), "NOT 0 <= deg(q["+str(j)+"]) < deg(q["+str(j-1)+"])"
        b[j] = b[j-2] + f[j]*b[j-1]
        if debug:
            print("\nИТЕРАЦИЯ НОМЕР "+str(j))
            print("Длелим q[j-2] = q["+str(j-2)+"] = "+str(q[j-2]))
            print("на q[j-1] = q["+str(j-1)+"] = "+str(q[j-1]))
            print("Получаем q[j-2] = q["+str(j-2)+"] = f[j]*q[j-1] + q[j] = f["+str(j)+"]*q["+str(j-1)+"] + q["+str(j)+"]")
            print("q["+str(j-2)+"] = ("+str(f[j])+")*("+str(q[j-1])+") + ("+str(q[j])+")")
            print("Находим b[j] = b["+str(j)+"] = b[j-2] + f[j]*b[j-1] = b["+str(j-2)+"] + f["+str(j)+"]*b["+str(j-1)+"]")
            print("b["+str(j)+"] = ("+str(b[j-2])+") + ("+str(f[j])+")*("+str(b[j-1])+") = "+str(b[j]))
    if debug:
        print("Так как deg(q[j]) = "+str(q[j].deg)+" <= t = "+str(t)+", то заканчиваем алгоритм")
    if type(b[j][0]) == GFEps:
        j += 1
        data["j"] = j
        b[j] = (b[j-1] / b[j-1][0])[0]
        if debug:
            print("Приводим b[j] = b["+str(j-1)+"] к виду полинома локаторов ошибок (должна быть единицу при 0 степени)")
            print("То есть, делим b["+str(j-1)+"] на "+str(b[j-1][0])+")")
            print("b["+str(j-1)+"]/"+str(b[j-1][0])+" = "+str(b[j]))
    data["Lambda"] = b[j]
    if return_data:
        return (data, b[j])
    return b[j]

def EuA2(S, t, debug = False):
    S_x = DPoly([1]+[S[i] for i in S])
    r = {0: DPoly(2*t+1), 1: S_x}
    b = {0: 0, 1: 1}
    a = {0: 1, 1: 0}
    q = {}
    j = 1
    if debug:
        print("Начальные значения:")
        print("S(x) =", S_x)
        print("r0 =", r[0])
        print("r1 =", r[1])
        print("b0 =", b[0])
        print("b1 =", b[1])
        print("a0 =", a[0])
        print("a1 =", a[1])
        print("j  =", j)
    while r[j].deg > t:
        j += 1
        q[j], r[j] = r[j-2]/r[j-1]
        a[j] = a[j-2] + q[j]*a[j-1]
        b[j] = b[j-2] + q[j]*b[j-1]
        if debug:
            print("Итерация номер", j)
            print("Делим r"+str(j-2)+" на r"+str(j-1))
            print("r"+str(j-2)+" = q"+str(j)+"*r"+str(j-1)+" + r"+str(j))
            print("r"+str(j-2)+" = ("+str(q[j])+")*("+str(r[j-1])+") + ("+str(r[j]))
            print("r"+str(j)+" = "+str(r[j]))
            print("q"+str(j)+" = "+str(q[j]))
            print("Находим a"+str(j)+" = a"+str(j-2)+" + q"+str(j)+"*a"+str(j-1))
            print("a"+str(j)+" = "+str(a[j]))
            print("Находим b"+str(j)+" = b"+str(j-2)+" + q"+str(j)+"*b"+str(j-1))
            print("b"+str(j)+" = "+str(b[j]))
    Omega  = r[j]
    Lambda = b[j]
    if debug:
        print("Так как deg(r"+str(j)+") <= t = "+str(t)+" - заканчиваем алгоритм")
        print("Omega  = r"+str(j)+" = "+str(r[j]))
        print("Lambda = b"+str(j)+" = "+str(b[j]))
    if Omega[0] != 1:
        Omega = (Omega/Omega[0])[0]
    if Lambda[0] != 1:
        Lambda = (Lambda/Lambda[0])[0]
    if debug:
        print("Приводим многочлены в необходимому виду:")
        print("Lambda = "+str(Lambda))
        print("Omega  = "+str(Omega))
    return (Lambda, Omega)

#Chen algorithm to find roots of error locator polynom
def ChenA(Lambda, field, debug = False):
    roots = []
    for i in range(0, len(field.field)):
        res = Lambda(field(-i))
        if debug:
           print("Возьмём степень ", i)
           print("Найдем Lambda(eps^"+str(-i)+")")
           print("Lambda(eps^"+str(-i)+") = "+str(res))
        if res == 0:
            if debug:
                print("Так как Lambda(eps^"+str(-i)+") = 0, то получаем один из корней полинома")
            roots.append(field(i))
    if debug:
        print("В итоге получаем корни:")
        print(roots)
    return roots

#Berlekamp-Massey algorithm for key equation
def BMA(S, t, eps_0, debug = False, return_data = False):
    try:
        r = 0
        Lambda = {0: eps_0}
        L      = {0: 0}
        B      = {0: DPoly([1])}
        Omega  = {0: DPoly([0])}
        A      = {0: DPoly([1])}
        dr     = {}
        M      = {}
        data   = {"r" : r, "S" : S, "dr" : dr, "M" : M, "B": B, "Lambda" : Lambda, "L" : L, "Omega" : Omega, "A" : A}
        while True:
            if debug:
                print("\nВывод значений для записи в таблицу")
                for d in data:
                    if d == "r":
                        print("r", "\t\t:\t", str(r))
                    else:
                        try:
                            print(d, "\t\t:\t", str(data[d][r]))
                        except KeyError:
                            print(d, "\t\t:\t---")
                print("\n\n")
            r += 1
            if debug:
                print("Увеличиваем r на еденицу")
                print("\nИТЕРАЦИЯ НОМЕР "+str(r))
            data["r"] = r
            if debug:
                print("Расчитываем dr по формуле: dr = Lambda_0 * S_r + ... + Lambda_L * S_r-L, то есть")
                print("Сумма от j=0 до L: Lambda_j*S_r-j, где Lambda_j - j-ый коэффициент многочлена Lambda(x)")
                print("S_r-j - r-j-ый синдром")
            dr[r] = Lambda[r-1][0]*S[r]
            for j in range(1, L[r-1]+1):
                dr[r] = dr[r] + Lambda[r-1][j]*S[r-j]
            if debug:
                print("dr = ", dr[r])
            if dr[r] == 0:
                if debug:
                    print("Так как dr = 0, то есть существуюший регистр генерирует следующий компонент синдрома")
                    print("то отводы обратной связи остаются без изменений:")
                    if Lambda:
                        print("Lambda(x) = " + str(Lambda[r-1]))
                    if M:
                        print("M(x)      = " + str(M[r-1]))
                    if L:
                        print("L         = " + str(L[r-1]))
                    print("И рассчитывается B(x) = x * B(x) = " + str(DPoly([0,1])) + " * " + str(B[r-1]))
                    
                B[r] = DPoly([0,1])*B[r-1]
                if debug:
                    print("B(x) = " + str(B[r]))
                if Lambda:
                    Lambda[r] = Lambda[r-1]
                if M:
                    M[r] = M[r-1]
                if L:
                    L[r] = L[r-1]
            else:
                if debug:
                    print("Так как dr != 0, то есть существующий регистр не генерирует следующий компонент синдрома")
                    print("то отводы обратной связи необходимо изменить:")
                    print("вычистить новый многочлен связей M(x) = Lambda(x) + dr*x*B(x)")
                s = dr[r]*(DPoly([0,1])*B[r-1])
                M[r] = s.__add__(Lambda[r-1], True)
                if debug:
                    print("M(x) = "+str(Lambda[r-1])+" + "+str(dr[r])+" * x * "+str(B[r-1]))
                    print("M(x) = "+str(M[r]))
                if 2*L[r-1] <= r-1:
                    if debug:
                        print("Так как 2*L <= r-1, то есть необходимо увеличить длину регистра, то")
                    B[r] = (dr[r]**-1)*Lambda[r-1]
                    Lambda[r] = M[r]
                    L[r] = r - L[r-1]
                    if debug:
                        print("рассчитываем B(x) = dr^-1 * Lambda(x) = "+str(dr[r]**-1)+" * "+str(Lambda[r-1])+" = "+str(B[r]))
                        print("Lambda(x) = M(x) = "+str(M[r]))
                        print("L = r-L = "+str(r)+" - "+str(L[r-1])+" = "+str(L[r]))
                else:
                    if debug:
                        print("Так как 2*L > r-1, то есть увеличивать длину регистра не надо, то")
                    L[r] = L[r-1]
                    Lambda[r] = M[r]
                    B[r] = DPoly([0,1])*B[r-1]
                    if debug:
                        print("L оставляем той же L = "+str(L[r]))
                        print("Lambda(x) = M(x) = "+str(Lambda[r]))
                        print("B(x) = x*B(x) = x * "+str(B[r-1])+" = "+str(B[r]))                    
            Omega[r] = Omega[r-1] + dr[r]*A[r-1]
            if debug:
                print("Расчитываем многочлен значений ошибок:")
                print("Omega_r(x) = Omega_r-1(x) + dr * A_r-1(x) = "+str(Omega[r-1])+" + "+str(dr[r])+" * "+str(A[r-1])+" = "+str(Omega[r]))
            if L[r] - L[r-1] == 0:
                A[r] = DPoly([0,1])*A[r-1]
                if debug:
                    print("Так как L_r - L_r-1 = 0, то расчитываем А(x):")
                    print("A(x) = x*A(x) = "+str(DPoly([0,1]))+" * "+str(A[r-1])+" = "+str(A[r]))
            else:
                A[r] = (dr[r]**-1)*DPoly([0,1])*Omega[r-1]
                if debug:
                    print("Так как L_r - L_r-1 != 0, то рассчитываем А(x):")
                    print("A(x) = dr^-1 * x * Omega_r-1(x) = " + str(dr[r]**-1)+" * x * "+str(Omega[r-1])+" = "+str(A[r]))
            if r == 2*t:
                if type(Lambda[r]) == DPoly:
                    if L[r] != Lambda[r].get_deg():
                        if debug:
                            print("В принятой комбинации больше t ошибок или ошибок нет вообще!")
                        return DPoly([1])
                else:
                    if L[r] != 0:
                        if debug:
                            print("В принятой комбинации больше t ошибок или ошибок нет вообще!")
                        return DPoly([1])
                if debug:
                    print("Так как r = 2*t = "+str(2*t)+", то алгоримт закончен и L = deg(Lambda(x)) = "+str(L[r]))
                break
            if debug:
                print("Так как r != 2*t продолжаем итерацию")
        if return_data:
            return (data, data["Lambda"][len(data["Lambda"])-1])
        return data["Lambda"][len(data["Lambda"])-1]
    except Exception as e:
        if debug:
            print("Возникла ошибка: "+ e)
            print("В принятой комбинации больше t ошибок или ошибок нет вообще!")
        return DPoly([1])

def BMA2(S, t, debug = True):
    Lambda = DPoly([1])
    ro     = DPoly([0,1])
    i      = 0
    l      = 0
    if debug:
        print("Начальные значения:")
        print("Lambda =", Lambda)
        print("ro =", ro)
        print("i =", i)
        print("l =", l)
    while i < 2*t + 1:
        i += 1
        if i >= 2*t + 1:
            break
        if debug:
            print("Увеличиваем i на единицу => i =", i)
        dr = Lambda[0]*S[i]            
        for j in range(1, l+1):
            dr += Lambda[j]*S[i-j]
        if debug:
            print("Ищем dr как сумма от j = 0 до l - Lambda[j]S[1-j]")
            print("dr =", dr)
        if dr != 0:
            Lambda_new = Lambda + dr*ro
            if debug:
                print("Так как dr != 0, то ищем Lambda_new = Lambda + dr*ro")
                print("Lambda_new =", Lambda_new)
            if 2*l < i:
                l = i - l
                ro = (Lambda/dr)[0]
                if debug:
                    print("Так как 2*l < i, то ищем l = i - l и ro = (Lambda/dr)")
                    print("l =", l)
                    print("ro =", ro)
            Lambda = Lambda_new
            if debug:
                print("Lambda = Lambda_new =", Lambda_new)
        ro = DPoly([0,1])*ro
        if debug:
            print("Находим ro = x * ro")
            print("ro =", ro)
    if debug:
        print("Так как i < 2*t - 1, то заканчиваем алгоритм")
    return Lambda


def find_G(powers, field):
    G = 1 + DPoly(1)*field(powers[0])
    for i in range(1, len(powers)):
        G *= 1 + DPoly(1)*field(powers[i])
    return G

def BMARS(G, S, t, debug = False):
    try:
        l = G.deg
        r = L = l
        Psi = B = G
        while r != 2*t:
            r += 1
            dr = S[r]
            summ = Psi[1]*S[r-1]
            for j in range(2, L+1):
                summ += Psi[j]*S[r-j]
            dr += summ
            if dr == 0:
                B = DPoly(1)*B
            else:
                M = Psi + dr*DPoly(1)*B
                if 2*L <= G.deg + r - 1:
                    B = (dr**-1)*Psi
                    Psi = M
                    L = r + G.deg - L
                else:
                    Psi = M
                    B = DPoly(1)*B
            if debug:
                print("\nr =",r)
                print("S[r] =", S[r])
                print("dr =", dr)
                print("L =",L)
                print("Psi =", Psi)
                print("B =", B)
                if M:
                    print("M =", M)
    except Exception as e:
        if debug:
            print(e)
            print("Ошибок нет или их больше t")
        return DPoly([1])
    return Psi

def find_syndromes(poly, field, l):
    S = {}
    for i in range(1, l+1):
        S[i] = poly(field(i))
        print(i)
        print(S[i])
    print(S)
    return S

def ForniA(roots, Omega, Lambda, t, debug = False):
    dLambda = Lambda.der()
    if debug:
        print("Находим производную Lambda(x)")
        print("Lambda'(x) =", dLambda)
    errors = []
    for root in roots:
        value = Omega(root**-1)*((dLambda(root**-1))**-1)
        errors.append(value)
        if debug:
            print("Рассчитываем Omega(Xi^-1)/Lambda'(Xi^-1):")
            print("Omega("+str(root**-1)+")/Lambda'("+str(root**-1)+") = "+str(value))
    return errors

def ForniA2(roots, Omega, Lambda, t, b, debug = False):
    dLambda = Lambda.der()
    if debug:
        print("Находим производную Lambda(x)")
        print("Lambda'(x) =", dLambda)
    errors = []
    for root in roots:
        value = (root**(1-b))*Omega(root**-1)*((dLambda(root**-1))**-1)
        errors.append(value)
        if debug:
            print("Рассчитываем Xi**(1-b) * Omega(Xi^-1)/Lambda'(Xi^-1):")
            print(Omega(root**-1))
            print(dLambda(root**-1))
            print(str(root**(1-b))+" * Omega("+str(root**-1)+")/Lambda'("+str(root**-1)+") = "+str(value))
    return errors


if __name__=="__main__":
    gf = GF([1,1,0,0,1])
    e = DPoly([0,gf(0), gf(1), gf(7), 0,0,0,0,0,0,gf(11)])
    S = {1: gf(0), 2: gf(13), 3: gf(12), 4: gf(5), 5: gf(0), 6: gf(3)}
    G = find_G([1,2], gf)
    Psi = BMARS(G, S, 3, True)
    S_x = DPoly([S[i] for i in S])
    Omega = ((S_x*Psi)/DPoly(6))[1]
    print(Psi)
    print(ForniA(ChenA(Psi, gf),Omega, Psi, 3, 1)) 
