import matplotlib.pyplot as pyplot
import numpy
import control.matlab as matlab
from numpy.lib.nanfunctions import nanargmax
import math
import colorama as color
from scipy.integrate import quad


w1 = matlab.tf([1], [4, 1])
w2 = matlab.tf([3], [6, 1])
w3 = matlab.tf([21], [4, 1])
# PID调节器：k=0.131, Ti=8, Tg=7; PI 调节器： k=0.03, Ti=13
k = 0.03
Ti = 13
# Tg = 7
# w_pid = matlab.tf([Ti*Tg*k, k*Ti, k], [Ti, 0])
w_pi = matlab.tf([k*Ti, k], [Ti, 0])
w_fenzi = w1*w2*w3*w_pi #对于PI调节器需要改这里
w_fenmu = 1+w1*w2*w3*w_pi
W = matlab.tf(w_fenzi.num, w_fenmu.num)
print(W)
# 阶跃响应
timeLine = []
for i in range(0, 10000): #对于ПИ的时候需要取到10000, PID 的时候跳到40000
    timeLine.append(i / 100)
[y1, x1] = matlab.step(W, timeLine)
pyplot.plot(x1, y1, 'green')
pyplot.ylabel('Амплитуда')
pyplot.xlabel('Время [c]')
pyplot.title('Переходная характеристика')
pyplot.grid()
pyplot.show()
long_y1 = int(len(y1))
y1_last = y1[int(len(y1)) - 1]
# print(long_y1)
# Интеграл
m = 0
M = 0
for i in range(long_y1 - 1):
    m = abs((y1[i]+y1[i+1])/2-y1_last)
    M += m * (x1[i+1]-x1[i])
# print(x1[3200]-x1[3199])

y1_max1 = y1[0]
y1_max2 = y1[0]
i_max1 = 0
i_max2 = 0

flag2 = True
list_4_0 = []
tp_list = []
flag1 = True
tp = 0
minx = 5
ymin = y1[0]
for i in range(long_y1):
    if i + 1 < long_y1:
        q_zuo = y1[i] - y1_last
        q_you = y1[i + 1] - y1_last
        if q_zuo * q_you <= 0 and flag1:
            list_4_0.append(i / 100)
            if len(list_4_0) > 9:
                flag1 = False
    if abs(y1[i]) > y1_max1:
        y1_max1 = y1[i]
        i_max1 = i

    if i > i_max1 and (i + 1) < long_y1 and i > 0:
        if y1[i] > y1[i + 1] and y1[i] > y1[i - 1] and flag2:
            y1_max2 = y1[i]
            i_max2 = i
            flag2 = False
    if abs(y1[i] - y1_last) >= 0.05 * y1_last:
        if i / 100 > tp:
            tp = i / 100
    # if abs(2.782608695652-i/100)<minx:
    #     minx = abs(2.782608695652-i/100)
    #     ymin = y1[i]

# print(color.Fore.BLUE +'ymin',y1[278])
print(list_4_0)
print(color.Fore.YELLOW + 'y1_last', y1_last)
# tk = list_4_0[2]
tk = (i_max2 - i_max1)/100
print('y1_max1', y1_max1)
print('t1', i_max1 / 100)
print('y1_max2', y1_max2)
print('t2', i_max2 / 100)
print('время регулирования:', tp)
print('колебательности:', tp / tk)
print('перерегулирование:', (y1_max1 - y1_last) / y1_last)
print('степень затухания:', (y1_max1 - y1_max2) / (y1_max1 - y1_last))

# распределению корней и полюсов на комплексной плоскости

poles_1, zeros = matlab.pzmap(W, plot=True, grid=True, title='Pole Zero Map')
pyplot.show()
poles = numpy.roots(W.den[0][0])
print(poles)

min_real = abs(poles[0].real)
max_fi = math.atan(abs(poles[0].imag) / abs(poles[0].real))
for i in range(len(poles)):
    if abs(poles[i].real) < min_real:
        min_real = abs(poles[i].real)
        a = poles[i]

    if (math.atan(abs(poles[i].imag) / abs(poles[i].real))) > max_fi:
        max_fi = math.atan(abs(poles[i].imag) / abs(poles[i].real))
        μ_max = abs(poles[i].imag) / abs(poles[i].real)

    if poles[i].real > 0:
        print('Система неустойчива')
        flag = False
        break

μ = math.tan(max_fi)
print(color.Fore.GREEN + 'max_fi:', max_fi)
print(color.Fore.GREEN + 'μ:', μ)
print('время регулирования:', 3.0 / min_real)
print('степень колебательности:', μ)
print('колебательности:', 1 / μ)
print('перерегулирование:', math.e ** (-math.pi / μ))
print('степень затухания:', (1 - math.e ** (-2 * math.pi / μ)))

# Bode график

y3_, y4_, w_bode = matlab.bode(W)
pyplot.show()
y3 = []
y4 = []
for i in range(int(len(y3_))):
    # для получения децибел(dB)
    y3.append(math.log10(y3_[i]) * 20)
    y4.append(y4_[i] / math.pi * 180)
    w_bode[i] = math.log10(w_bode[i])
a = y3_[0]
b = y4_[0]
i_indekc1 = 0
i_indekc2 = 0
for i in range(int(len(y3_))-1):
    if y3[i] * y3[i+1] <= 0:
        i_indekc1 = i
for i in range(int(len(y4_))):
    if abs(y4[i]+180) <= 1: # 这里的精度对于PID等于4可以变，这里取4只是方便此题计算
        i_indekc2 = i

Aw = y3_[i_indekc2]
Kg = -1 * y3[i_indekc2]
Y = y4[i_indekc1] + 180
print(color.Fore.CYAN + 'Aw:', Aw)
print(color.Fore.CYAN + 'Kg/dB:', Kg)
print(color.Fore.CYAN + 'γ:', Y)


# АЧХ
def годограф(W):
    '''
    для пострения годографа
    '''
    wLine = []
    re = []
    im = []
    amplitude = []
    for i in range(0, 500):
        a = 0
        b = 0
        w = i / 100
        wLine.append(w)
        onepatop = w * 1j
        for i1 in range(len(W.num[0][0])):
            a += W.num[0][0][i1] * onepatop ** (int(len(W.num[0][0])) - i1 - 1)

        for i1 in range(len(W.den[0][0])):
            b += W.den[0][0][i1] * onepatop ** (int(len(W.den[0][0])) - i1 - 1)
        c = a / b
        re.append(c.real)
        im.append(c.imag)
        amplitude.append(math.sqrt(c.real ** 2 + c.imag ** 2))
    pyplot.plot(wLine, amplitude)
    pyplot.ylabel('A(w)')
    pyplot.xlabel('w[rad]')
    pyplot.title('АЧХ')
    pyplot.grid(True)
    pyplot.show()
    return wLine, re, im, amplitude

wLine, re, im, amplitude = годограф(W)

A_0 = amplitude[0]
A_max = amplitude[0]
wucha = 20
i_max = 0
wcp = 0
for i in range(int(len(wLine))):
    if amplitude[i] > A_max:
        A_max = amplitude[i]
        i_max = i
    if i > i_max and abs(amplitude[i] - A_0) < wucha:
        wcp = wLine[i]
        wucha = abs(amplitude[i] - A_0)

print(color.Fore.RED + 'wcp:', wcp)
print(color.Fore.RED + 'A_max:', A_max)
print(color.Fore.RED + 'A_0:', A_0)
print('Показатель колебательности:', A_max / A_0)
print('время регулирования:', 2*math.pi / wcp, '-', 2 *2* math.pi / wcp)

# Интеграл
print(color.Fore.BLUE + 'Интеграл:', M)
