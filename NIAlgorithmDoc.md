# Error Control codes algorithms
## There are realized four algorithms:
### Berlekamp–Massey algorithm | Алгоритм Берлекэмпа-Мэсси для решения ключевого уравнения
+ [Wiki](https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm)
+ [Вики](https://ru.wikipedia.org/wiki/%D0%90%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC_%D0%91%D0%B5%D1%80%D0%BB%D0%B5%D0%BA%D1%8D%D0%BC%D0%BF%D0%B0_%E2%80%94_%D0%9C%D1%8D%D1%81%D1%81%D0%B8)

```
BMA(S, t, eps_0, debug = False)
```
+ S - the dictionary of syndromes started with the index of 1: {1: S1, 2: S2, ..., n: Sn}
+ t - the minimum amount of error the code can detect
+ eps_0 - the element of Galois field (see [Galois Field realization](GFDoc.md)) of power 0 - eps^0
+ debug - the flag to see full information about the algorithm proccess (in russian language)

### Piterson-Gorenshtein-Zirler algorith | Алгоритм Питерсона-Горенштейна-Цирлера для решения ключевого уравнения
```
PGCA(S, t, degub = False)
```
+ S - the dictionary of syndromes started with the index of 1: {1: S1, 2: S2, ..., n: Sn}
+ t - the minimum amount of error the code can detect
+ debug - the flag to see full information about the algorithm proccess (in russian language)

### Euclid algorith to solve key equation | Алгоритм Евклида для решения ключевого уравнения
+ [Вики](https://ru.wikipedia.org/wiki/%D0%90%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC_%D0%95%D0%B2%D0%BA%D0%BB%D0%B8%D0%B4%D0%B0)
+ [Wiki](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm)
```
EuA(S, t, debug = False)
```
+ S - the dictionary of syndromes started with the index of 1: {1: S1, 2: S2, ..., n: Sn}
+ t - the minimum amount of error the code can detect
+ debug - the flag to see full information about the algorithm proccess (in russian language)

### Chen algorith to find roots of error locator polynomial | Алгоритм Ченя для нахождения корней полинома локаторов ошибок
```
ChenA(Lambda, field, debug = False)
```
+ Lambda - error locator polynomial of type DPoly (see [GF documentations](GFDoc.md))
+ field  - Galois field (also see [GF documentations](GFDoc.md))
+ debug - the flag to see full information about the algorithm proccess (in russian language)
