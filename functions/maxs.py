import numpy as np


def find_max(x, y, th=10, add=False, number:int=2, **kwargs):
    if th!=10:
        th -= 0.4 if not add else -0.4
    m = min(y)
    maxs = []
    for idx, i in enumerate(x):
        back, fowr = idx-1, idx+1
        back_2, fowr_2 = idx-2, idx+2
        if back<0: back_2=-2; back=-1
        if back_2<0 and back>=0: back_2=-1
        if fowr>=len(x): fowr=0; fowr_2=1
        if fowr_2==len(x): fowr_2=0 
        if y[idx] > y[back] and y[idx] > y[fowr] and y[idx] > m+th and y[idx] > y[back_2] and y[idx] > y[fowr_2] and (np.abs(y[idx]-y[fowr])>th/627.51 or np.abs(y[idx]-y[back])>th/627.51): maxs.append(i)
    if len(maxs) != number and 0.4 < th < 25:
        flag = len(maxs) < number
        if th == 10: th -= 0.4 if flag else -0.4
        maxs = find_max(x, y, th, add=not flag)
    return maxs




def find_max(x, y, th:float=10, add:bool=False, number:int=2, deep:int=4, hartree=False, red:float=0.2):



    # reduce energy treshold 
    if th!=10:
        th -= red if not add else -red
    # if th<0.2:
    #     print(f'Not found {number} max(s). Going to exit')
    #     return None

    m = min(y) +  th / (1 if not hartree else 627.51)
    maxs = []
    for idx, i in enumerate(list(y)):

        flag = False

        if (idx-deep)%len(x) < (idx+deep)%len(x):
            for j in range((idx-deep)%len(x), (idx+deep)%len(x)):
                flag = i > y[j] and i >= m
        else: 
            rng = list(range((idx-deep)%len(x), len(x)-1)) + list(range(0, (idx+deep)%len(x)))
            for j in rng:
                flag = i > y[j] and i >= m

        if flag: maxs.append(x[idx])


    if len(maxs) != number and (0.2 <= th):
        f = len(maxs) < number
        if th == 10: th -= (red if f else -red)
        maxs = find_max(x, y, th, add=not f, hartree=hartree)
    
    return maxs



if __name__ == '__main__':
    with open('tests/clockwise.log') as f:
        fl = f.readlines()
    x = [float(i.split()[0]) for i in fl]
    y = [float(i.split()[1]) for i in fl]
    m = find_max(x, y, th=10, hartree=True, number=2, deep=2)
    assert m == [215.2826, 400.42545714], f"NO. Obtained: {m}"