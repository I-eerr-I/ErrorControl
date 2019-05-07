class DMatrix:

    def __init__(self, matrix):
        self.matrix = [[matrix[i][j] for j in range(len(matrix[0]))]
                       for i in range(len(matrix))]
        self.update()
        
    def __getitem__(self, line):
        return self.matrix[line]

    def __setitem__(self, line, item):
        self.matrix[line] = item

    def __delitem__(self, item):
        del self.matrix[item]
        self.update()

    def update(self):
        self.tp = type(self.matrix[0][0])
        self.rows = len(self.matrix)
        self.columns = len(self.matrix[0])
        self.shape = (self.rows, self.columns)

    def copy(self):
        return DMatrix(self.matrix)

    def append(self, item):
        tp = self._check_type(item)
        if tp >= 0:
            self._check_shapes(item)
            if tp > 0:
                if not type(item[0]) == list:
                    for i in range(self.rows):
                        self[i].append(item[i])
                    return True
            for i in range(self.rows):
                for j in range(self.columns):
                    self[i].append(item[i][j])
            return True
        else:
            for i in range(self.rows):
                self[i].append(item)
            return True
                
    def get_zero(self, rows=None, columns=None):
        if not rows:
            rows = self.rows
        if not columns:
            columns = self.columns
        return DMatrix([[0 for _ in range(self.columns)]
                       for _ in range(self.rows)])

    def shape(self):
        return self.shape

    def diag(self):
        return [self.matrix[i][i] for i in range(self.rows)]

    def row(self, j):
        return [self.matrix[i][j] for i in range(self.rows)]

    def trans(self):
        return DMatrix([self.row(i) for i in range(self.columns)])

    def mminor(self, line, row=0):
        result = DMatrix(self.matrix)
        del result[line]
        for i in range(result.rows):
            for j in range(result.columns):
                if j == row:
                    c = result[i].copy()
                    del c[j]
                    result[i] = c
        result.update()
        return result

    def alg_adj(self, i, j=0):
        return self.mminor(i,j).det()*((-1)**(i+j))

    def malg_adj(self):
        result = DMatrix(self.matrix)
        for i in range(self.rows):
            for j in range(self.columns):
                result[i][j] = self.alg_adj(i, j)
        return result

    def inverse(self):
        self.update()
        if self.shape == (1,1):
            return DMatrix([self[0][0]**-1])
        d = self.det()
        if d == 0:
            return False
        result = DMatrix(self.matrix)
        algT = self.malg_adj().trans()
        for i in range(self.rows):
            for j in range(self.columns):
                result[i][j] = algT[i][j]*(abs(d)**-1)
        return result

    def det(self):
        if not self.shape[0] == self.shape[1]:
            raise Exception("DMatrix need to be squared")
        result = 0
        if self.shape == (1,1):
            return self.matrix[0][0]
        for i in range(self.rows):
            result += self.matrix[i][0]*self.alg_adj(i)
        return result
    
    def __list__(self):
        return [i for i in range(self.rows) for j in range(self.columns)]

    def __len__(self):
        return sum([len(self.matrix[i]) for i in range(self.rows)]) 

    def __neg__(self):
        return DMatrix([[-self[i][j] for j in range(self.columns)]
                       for i in range(self.rows)])
    
    def _check_shapes(self, other):
        if self._check_type(other) > 0:
            if (not self.rows == len(other)
                and self.columns == len(other[0])):
                raise Exception("Not equal shapes")
        elif not self.shape == other.shape:
            raise Exception("Not equal shapes")

    def _check_mul_shapes(self, other):
        if self._check_type(other) > 0:
            if not self.columns == len(other):
                raise Exception("Not right shapes")
        elif not self.columns == other.rows:
            raise Exception("Not right shapes")

    def _check_type(self, other):
        if type(other) == type(self):
            return 0
        elif type(other) == list:
            return 1
        else: return -1
    
    def __add__(self, other):
        tp = self._check_type(other)
        if tp >= 0:
            self._check_shapes(other)
            if tp > 0:
                other = DMatrix(other)
            return DMatrix([[self[i][j]+other[i][j] for j in range(self.columns)]
                  for i in range(self.rows)])
        else:
            return DMatrix([[self[i][j]+other for j in range(self.columns)]
                           for i in range(self.rows)])
        
    def __mul__(self, other):
        tp = self._check_type(other)
        if tp >= 0:
            self._check_mul_shapes(other)
            if tp > 0:
                other = DMatrix(other)
            return DMatrix([[sum([self[i][k]*other[k][j] for k in range(self.columns)])
                             for j in range(other.columns)]
                            for i in range(self.rows)])
        else:
            return DMatrix([[self[i][j]*other for j in range(self.columns)]
                           for i in range(self.rows)])
        
    def __pow__(self, power):
        result = DMatrix(self.matrix)
        if power == -1:
            return self.inverse()
        for k in range(power-1):
            result = result * self
        return result

    def __sub__(self, other):
        tp = self._check_type(other)
        if tp >= 0:
            self._check_shapes(other)
            if tp > 0:
                other = DMatrix(other)
            return DMatrix([[self[i][j]-other[i][j] for j in range(self.columns)]
                  for i in range(self.rows)])
        else:
            return DMatrix([[self[i][j]-other for j in range(self.columns)]
                           for i in range(self.rows)])
    def __str__(self):
        M_str = "" 
        for line in self.matrix:
            M_str += str(line) + "\n"
        return M_str.strip()

    def __repr__(self):
        return self.__str__()

    def __and__(self, other):
        tp = self._check_type(other)
        if tp >= 0:
            self._check_shapes(other)
            if tp > 0:
                other = DMatrix(other)
            return DMatrix([[self[i][j]&other[i][j] for j in range(self.columns)]
                         for i in range(self.rows)])
        else:
            return DMatrix([[self[i][j]&other for j in range(self.columns)]
                           for i in range(self.rows)])

    def __or__(self, other):
        self._check_shapes(other)
        return DMatrix([[self[i][j]|other[i][j] for j in range(self.columns)]
                       for i in range(self.rows)])

    def __eq__(self, other):
        return self.matrix == other.matrix

    def sign(self):
        return DMatrix([[1 if self[i][j] > 0 else 0
                        for j in range(self.columns)]
                           for i in range(self.rows)])

    def boolmul(self, other, times=0):
        from functools import reduce
        self._check_mul_shapes(other)
        result = DMatrix([[reduce(lambda x,y: x or y,[self[i][k] & other[k][j] for k in range(other.columns)])
                             for j in range(other.columns)]
                            for i in range(self.rows)])
        if times > 0:
            for k in range(times-1):
                result = DMatrix([[reduce(lambda x,y: x or y,[result[i][k] & other[k][j] for k in range(other.columns)])
                                 for j in range(other.columns)]
                                for i in range(result.rows)])
        return result



if __name__ == '__main__':
    test_suit = (
        [[5,6],
         [8,9]],
        [[4,6],
         [7,9]],
        [[4,5],
         [7,8]],
        [[2,3],
         [8,9]],
        [[1,3],
         [7,9]],
        [[1,2],
         [7,8]],
        [[2,3],
         [5,6]],
        [[1,3],
         [4,6]],
        [[1,2],
         [4,5]])
    testsuit_dets = (-3, -6, -3, -6, -12, -6, -3, -6, -3)
    A = DMatrix([[1,2,3],[4,5,6],[7,8,9]])
    k = 0
    for i in range(A.rows):
        for j in range(A.columns):
            minor = A.mminor(i,j)
            assert minor.matrix == test_suit[k], str(k)
            assert minor.det() == testsuit_dets[k], str(k) + ' ' + str(minor.det()) 
            k += 1 
