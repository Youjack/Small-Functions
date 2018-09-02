# -*- coding: utf-8 -*-

from NumListReader import ReadNumList

# Ext = [ Ext(x), Ext(y), Ext(z), ... ]

def ConvertToLinCoo(Ext, Coo):
    l = len(Ext)
    lExt = [ 1, ]
    for i in range(l-1):
        p = 1
        for j in range(i+1): p *= Ext[j]
        lExt.append(p)
    lCoo = 0
    for i in range(l):
        lCoo += Coo[i] * lExt[i]
    return lCoo

def ConvertFromLinCoo(Ext, lCoo):
    pass

if __name__ == '__main__':
    Ext = ReadNumList(input('Ext = '))
    Coo = ReadNumList(input('Coo = '))
    print(ConvertToLinCoo(Ext, Coo))
    input()
