# -*- coding: utf-8 -*-

# LinCoo < (Straight) Line Coordinate

from NumListReader import ReadNumList

# Ext = [ Ext(x), Ext(y), Ext(z), ... ]
# Coo = [ Coo(x), Coo(y), Coo(z), ... ]

def GenerateLinExt(d, Ext):
    lExt = [ 1, ]
    for i in range(d-1):
        p = 1
        for j in range(i+1): p *= Ext[j]
        lExt.append(p)
    return lExt

def ConvertToLinCoo(Ext, Coo : 'must be in Ext'):
    d = len(Ext)
    lExt = GenerateLinExt(d, Ext)
    lCoo = 0
    for i in range(d):
        lCoo += Coo[i] * lExt[i]
    return lCoo

def ConvertFromLinCoo(Ext, lCoo : 'must be in Ext'):
    d = len(Ext)
    lExt = GenerateLinExt(d, Ext)
    Coo = []
    for i in range(d): Coo.append(0)
    for i in range(d):
        Coo[d-i-1] = lCoo // lExt[d-i-1]
        lCoo %= lExt[d-i-1]
    return Coo

if __name__ == '__main__':
    toDo = int(input('To =0 or From =1 : '))
    if toDo == 0:
        Ext = ReadNumList(input('Ext = '))
        Coo = ReadNumList(input('Coo = '))
        print(ConvertToLinCoo(Ext, Coo))
    else:
        Ext = ReadNumList(input('Ext = '))
        lCoo = int(input('lCoo = '))
        print(ConvertFromLinCoo(Ext, lCoo))
    input()
