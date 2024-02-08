# -*- coding: utf-8 -*-

def ReadNumList(string):
    blank = ( ' ', ',', ';' )
    List = []
    s = ''
    for i in range(len(string)):
        if string[i] in blank:
            List.append(int(s))
            s = ''
        elif i == len(string)-1:
            List.append(int(s+string[i]))
        else:
            s += string[i]
    return List

if __name__ == "__main__":
    print(ReadNumList(input()))
    input()
