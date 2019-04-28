import matrix
from GF import DPoly, GF, GFEps

#Piterson-Gorenshtein-Zirler Algorithm for key equation
def PGCA(S, t, debug = False):
    Q   = 0
    det = 0
    Ms  = []
    for q in range(t, 0, -1):
        Ms = matrix.DMatrix([[S[i] for i in range(j, q+j)] for j in range(1, q+1)])
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
        return None
    Mi = Ms.inverse()
    M2 = matrix.DMatrix([[S[i] for i in range(q+1, 2*q+1)]])
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
    return data

#Euclid algorithm for key equation
def EuA(S, t, debug = False):
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
    return data

#Chen algorithm to find roots of error locator polynom
def ChenA(Lambda, field, debug = False):
    roots = []
    for i in range(0, len(field.field)):
        res = field.subs_value(Lambda, -i)
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
def BMA(S, t, eps_0, debug = False):
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
                    print(r, "\t\t:\t", str(r))
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
                print("Lambda(x) = " + str(Lambda[r-1]))
                print("M(x)      = " + str(M[r-1]))
                print("L         = " + str(L[r-1]))
                print("И рассчитывается B(x) = x * B(x) = " + str(DPoly([0,1])) + " * " + str(B[r-1]))
                
            B[r] = DPoly([0,1])*B[r-1]
            if debug:
                print("B(x) = " + str(B[r]))
            Lambda[r] = Lambda[r-1]
            M[r] = M[r-1]
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
            L = Lambda[r].get_deg()
            if debug:
                print("Так как r = 2*t = "+str(2*t)+", то алгоримт закончен и L = deg(Lambda(x)) = "+str(L))
            break
        if debug:
            print("Так как r != 2*t продолжаем итерацию")
            
    return data

if __name__=="__main__":
    S = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
    PGCA(S, 3, True)
    gf = GF(4, [1,1,0,0,1])
    ChenA(DPoly([1, gf(12), gf(13)]), gf, True)

