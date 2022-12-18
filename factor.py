#------------------------------------------------------------------------------
# # MAT0122 ÁLGEBRA LINEAR I
# Aluno: RODRIGO DE CASTRO MICHELASSI
# Numero USP: 13672703
# Tarefa: EP1 - Fatoração de Inteiros
# Data: 15/12/2022
# 
# DECLARO QUE SOU O ÚNICO AUTOR E RESPONSÁVEL POR ESTE PROGRAMA.  TODAS AS 
# PARTES DO PROGRAMA, EXCETO AS QUE SÃO BASEADAS EM MATERIAL FORNECIDO  
# PELO PROFESSOR OU COPIADAS DO LIVRO OU DO MATERIAL DIDÁTICO DE MAT0122, 
# FORAM DESENVOLVIDAS POR MIM.  DECLARO TAMBÉM QUE SOU RESPONSÁVEL POR TODAS 
# AS CÓPIAS DESTE PROGRAMA E QUE NÃO DISTRIBUÍ NEM FACILITEI A DISTRIBUIÇÃO
# DE CÓPIAS DESTA PROGRAMA.
#------------------------------------------------------------------------------

import sys
import GF2
from GF2 import one
from vec import Vec
from vecutil import list2vec
from factoring_support import intsqrt, dumb_factor, primes, prod, gcd
from functools import reduce    # used for the "find a and b" function
import operator

# returns a matrix on the echelon form and the amount Z of 0-vector rows in it
def transformation_rows(rowlist_input, col_label_list = None):
    z = 0
    one = GF2.one
    rowlist = list(rowlist_input)
    if col_label_list == None: 
        col_label_list = sorted(rowlist[0].D, key=repr)
    m = len(rowlist)
    row_labels = set(range(m))
    M_rowlist = [Vec(row_labels, {i:one}) for i in range(m)]
    new_M_rowlist = []
    rows_left = set(range(m))
    for c in col_label_list:
        rows_with_nonzero = [r for r in rows_left if rowlist[r][c] != 0]
        if rows_with_nonzero != []:
            pivot = rows_with_nonzero[0]
            rows_left.remove(pivot)
            new_M_rowlist.append(M_rowlist[pivot])
            for r in rows_with_nonzero[1:]:
                multiplier = rowlist[r][c] / rowlist[pivot][c]
                rowlist[r] -= multiplier * rowlist[pivot]
                M_rowlist[r] -= multiplier * M_rowlist[pivot]

    for r in rows_left: 
        new_M_rowlist.append(M_rowlist[r])
        z += 1
            
    return new_M_rowlist, z

def int2GF2(num):
    if num % 2 == 1:
        return GF2.one
    return 0

def prod(factors):
    return reduce(operator.mul, factors, 1)

def find_candidates(num, primeSet):
    factorList = []
    roots = []
    rowlist = []
    count = 0
    addX = 2
    
    while count < len(primeSet) + 1:
        a = {}
        x = intsqrt(num) + addX        
        factorList = dumb_factor(x * x - int(num), primeSet)
        if len(factorList) != 0:
            for i, j in factorList:
                value = j % 2 if j % 2 != 1 else GF2.one
                a[i] = value
            rowlist.append(Vec(primeSet, a))
            roots.append(x)
            count += 1

        addX += 1
        factorList.clear()                      # remove everything from factorList so it can be used again on next iteration

    # rowlist represents factorList vector over GF2 produced by make_vec
    return roots, rowlist

def find_a_and_b(v, roots, num):
    i = 0
    alist = []
    auxList = []

    while i < len(roots):
        if v[i] != 0:
            alist.append(roots[i])
        i += 1

    a = prod(alist)
    for x in alist:
        auxList.append(x**2 - num)

    c = prod(auxList)
    b = intsqrt(c)
    if b**2 == c:
        return (a,b)
    return (a, 0)

def findGCD(num, auxMatrix, paramList, roots):
    i = 0
    a = b = num
    gcdCalc = gcd(a - b, num)
    while (gcdCalc == 1 or gcdCalc == num) and i < len(auxMatrix):
        a, b = find_a_and_b(auxMatrix[i], roots, num)
        if len(paramList) >= 4 and b != 0:
            print(i, ": a = ", a, "/ b = ", b)
        if b != 0:
            gcdCalc = gcd(int(a) - int(b), num)
        i += 1
    return gcdCalc
    
def main():
    # getting params from the command line and defining U
    paramList = []
    rowList = []
    roots = []
    echelonMatrix = []
    auxMatrix = []
    primeSet = ()
    i = 0
    num = int(sys.argv[1])

    for param in sys.argv:
        paramList.append(param)

    if len(paramList) < 3:
        U = 10000
    else:
        U = int(paramList[2])
            
    # defining a set with the U first prime numbers
    primeSet = primes(U)

    roots, rowList = find_candidates(num, primeSet)
    echelonMatrix, z = transformation_rows(rowList)
    auxMatrix = echelonMatrix
    auxMatrix.reverse()

    gcdCalc = findGCD(num, auxMatrix, paramList, roots)

    if gcdCalc != 1 and gcdCalc != num:
        print("factor = ", gcdCalc)
    else:
        print("Failed")
    exit
main()