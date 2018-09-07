# -*- coding: utf-8 -*-

# Only Support Base-2 through Base-10

def ConvertBaseFrom(base : '2~10', num) -> int:
    num = str(num)
    _num = 0
    l = len(num)
    for i in range(l):
        _num += int(num[i]) * base ** (l-i-1)
    return _num

def ConvertBaseTo(base : '2~10', num) -> str:
    num = int(num)
    _num = ''
    while (num != 0):
        _num = str(num % base) + _num
        num = num // base
    return _num

if __name__ == '__main__':
    toDo = int(input('From =0 or To =1 ? : '))
    base = int(input('base = '))
    num = input('num = ')
    if toDo == 0: print(ConvertBaseFrom(base, num))
    else: print(ConvertBaseTo(base,num))
    input()
