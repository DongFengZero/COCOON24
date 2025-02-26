import math
import matplotlib.pyplot as plt
###input para(S>mn)
font = {'family': 'serif', 'serif': ['Times New Roman'], 'size': 15, 'style': 'normal'}
plt.rc('font', **font)
x_data=[]
y_data = []
y_data2 = []
fig, axes = plt.subplots(1,3)
max_temp=0
for i in range(20):
    m = 50 + 10 * i
    n = int(2500 / m)
    S = 2500
    x_data.append(m)

    R1S = min(m*n,S)
    R2S = min(m*n,S)

    T1S = m * m
    T2S = m*m
    T3S = m*m
    X = max(m * n,S)

    for k in range(n):
        Sre0 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
        T1S = min(math.floor(max(m * n - (k - 1) * (n - (k - 1) / 2) - S, m)), T1S)

        # Sre1 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
        Sre1 = Sre0 - (m - k + 1)
        Srem1 = Sre0 - (3 * m - 3 * k + 1)
        T21S = max(3 * (m - k) + 1 - max((math.ceil(Sre0)), 0), 0)
        T31S = max((m - k + 1) - max((math.ceil(Sre0)), 0), 0)

        Srem2 = Srem1 - (5 * m - 5 * k + 2) * (n - k + 1)
        Sre2 = Sre1 - (m - k + 1) * (n - k + 1)
        T22S = max((5 * m - 5 * k + 2) * (n - k + 1) - max(math.ceil(Srem1), 0), 0)
        T32S = max((m - k + 1) * (n - k + 1) - max(math.ceil(Sre1), 0), 0)

        Sre3 = Sre2
        Srem3 = Srem2
        T23S = max(2 * (m - k + 1) * (n - k + 1) - S, 0)
        T33S = max((m - k + 1) * (n - k + 1) - S, 0)

        T24S = max(m * n - k * (n - k / 2) - S, 0) + (m + n + 1) + max(2 * m * k - math.ceil(Srem3), 0)
        T34S = max(m * n - k * (n - k / 2) - S, 0)

        T2S = min(T21S + T22S + T23S + T24S, T2S)
        T3S = min(T31S + T32S + T33S + T34S, T3S)


    TS = 0
    if R1S-T2S <= R2S-T1S-T3S:
        TS = X-R2S+T1S+T3S
    else:
        TS = X-R1S+T2S

    TS1 = X - max(R1S,R2S) + min(T1S,T2S)

    V_1 = 4 * m * n + (4 * n - 4 * n * n) / 2
    V_2 = (3*m*n*n+5*m*n-2*n*n-n*n*n+n) / 2
    V_3 = (n - n * n * n) / 6 + (m * n + m * n * n) / 2
    V_4 = 5*n / 6 - m * n + m * n * n +n*n/2- n * n * n / 3
    V_5 = m * n * n - (n - 1) * n * (4 * n + 1) / 12
    V_all = math.ceil(V_1 + V_2 + V_3 + V_4 + V_5)

    h1 = m + 1 - max(math.floor(n + 1 - math.sqrt(max(n * n - 2 * m * n + 2 * S, 0))), 1)
    phi_1 = 4 * h1
    Phi_1 = h1

    h2 = S + Phi_1
    k2 = m + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
    i2 = n + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
    phi_2 = 3 * k2 * i2 - i2 + k2
    Phi_2 = k2 * i2

    h3 = S + Phi_1 + Phi_2
    k_3 = max(min(math.floor(
        (-math.sqrt(max(m * m + 4 * n * n + 1 + 4 * n - 8 * m * n - 4 * m + 6 * h3, 0)) + m + 2 * n + 4) / 3), n), 1)
    phi_3 = (m - k_3 + 1) * (n - k_3 + 1)
    Phi_3 = (m - k_3 + 1) * (n - k_3 + 1)

    h_4 = (k_3 - 1) * (k_3 - 1) / 2 + k_3 * (m - k_3 + 1) + (m - k_3) * (n - k_3)
    k_4 = max(min(math.floor(n - math.sqrt(max(n * n - n * m * 2 + 2 * h_4 - 1, 0))), n), 1)
    phi_4 = -k_4 * k_4 + (2 * m + 2) * k_4 - 2 * m

    V_max = S+phi_1+phi_2+phi_3+phi_4
    print(m,n,S,TS/TS1)
    y_data.append(TS1 * (V_all / V_max - 1))
    y_data2.append(TS * (V_all / V_max - 1))
    max_temp = max(TS / TS1, max_temp)

print(max_temp," ")
axes[0].grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
axes[0].grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)
axes[0].set_title("(a)S=2500,M>=N,N=INT(S/M)", fontsize=15,font='Times New Roman')
axes[0].plot(x_data,y_data2,color='red',linewidth=2.5,linestyle='--',label='(X1,X2)-Partition',marker='o')
axes[0].plot(x_data,y_data,color='blue',linewidth=2.0,linestyle='-.',label='X-Partition',marker='^')
axes[0].set_xlabel('The Size of M of Matrix', fontsize=15)
axes[0].set_ylabel('Q(I/O Complexity)', fontsize=15)
axes[0].legend(loc="best", fontsize=15)

y_data=[]
y_data2=[]
x_data=[]
max_temp=0
for i in range(45):
    m = 100 + 20 * i
    n = int(10000 / m)
    S = 10000
    x_data.append(m)

    R1S = min(m * n, S)
    R2S = min(m * n, S)

    T1S = m * m
    T2S = m * m
    T3S = m * m
    X = max(m * n, S)

    for k in range(n):
        Sre0 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
        T1S = min(math.floor(max(m * n - (k - 1) * (n - (k - 1) / 2) - S, m)), T1S)

        # Sre1 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
        Sre1 = Sre0 - (m - k + 1)
        Srem1 = Sre0 - (3 * m - 3 * k + 1)
        T21S = max(3 * (m - k) + 1 - max((math.ceil(Sre0)), 0), 0)
        T31S = max((m - k + 1) - max((math.ceil(Sre0)), 0), 0)

        Srem2 = Srem1 - (5 * m - 5 * k + 2) * (n - k + 1)
        Sre2 = Sre1 - (m - k + 1) * (n - k + 1)
        T22S = max((5 * m - 5 * k + 2) * (n - k + 1) - max(math.ceil(Srem1), 0), 0)
        T32S = max((m - k + 1) * (n - k + 1) - max(math.ceil(Sre1), 0), 0)

        Sre3 = Sre2
        Srem3 = Srem2
        T23S = max(2 * (m - k + 1) * (n - k + 1) - S, 0)
        T33S = max((m - k + 1) * (n - k + 1) - S, 0)

        T24S = max(m * n - k * (n - k / 2) - S, 0) + (m + n + 1) + max(2 * m * k - math.ceil(Srem3), 0)
        T34S = max(m * n - k * (n - k / 2) - S, 0)

        T2S = min(T21S + T22S + T23S + T24S, T2S)
        T3S = min(T31S + T32S + T33S + T34S, T3S)

    TS = 0
    if R1S - T2S <= R2S - T1S - T3S:
        TS = X - R2S + T1S + T3S
    else:
        TS = X - R1S + T2S

    TS1 = X - max(R1S, R2S) + min(T1S, T2S)

    V_1 = 4 * m * n + (4 * n - 4 * n * n) / 2
    V_2 = (3*m*n*n+5*m*n-2*n*n-n*n*n+n) / 2
    V_3 = (n - n * n * n) / 6 + (m * n + m * n * n) / 2
    V_4 = 5*n / 6 - m * n + m * n * n +n*n/2- n * n * n / 3
    V_5 = m * n * n - (n - 1) * n * (4 * n + 1) / 12
    V_all = math.ceil(V_1 + V_2 + V_3 + V_4 + V_5)

    h1 = m + 1 - max(math.floor(n + 1 - math.sqrt(max(n * n - 2 * m * n + 2 * S, 0))), 1)
    phi_1 = 4 * h1
    Phi_1 = h1

    h2 = S + Phi_1
    k2 = m + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
    i2 = n + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
    phi_2 = 3 * k2 * i2 - i2 + k2
    Phi_2 = k2 * i2

    h3 = S + Phi_1 + Phi_2
    k_3 = max(min(math.floor(
        (-math.sqrt(max(m * m + 4 * n * n + 1 + 4 * n - 8 * m * n - 4 * m + 6 * h3, 0)) + m + 2 * n + 4) / 3), n), 1)
    phi_3 = (m - k_3 + 1) * (n - k_3 + 1)
    Phi_3 = (m - k_3 + 1) * (n - k_3 + 1)

    h_4 = (k_3 - 1) * (k_3 - 1) / 2 + k_3 * (m - k_3 + 1) + (m - k_3) * (n - k_3)
    k_4 = max(min(math.floor(n - math.sqrt(max(n * n - n * m * 2 + 2 * h_4 - 1, 0))), n), 1)
    phi_4 = -k_4 * k_4 + (2 * m + 2) * k_4 - 2 * m
    V_max = S + phi_1 + phi_2 + phi_3 + phi_4
    print(m,n,S,TS/TS1)
    y_data.append(TS1 * (V_all / V_max - 1))
    y_data2.append(TS * (V_all / V_max - 1))
    max_temp=max(TS/TS1,max_temp)

print(max_temp," ")
axes[1].grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
axes[1].grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)
axes[1].set_title("(b)S=10000,M>=N,N=INT(S/M)", fontsize=15)
axes[1].plot(x_data,y_data2,color='red',linewidth=2.5,linestyle='--',label='(X1,X2)-Partition',marker='o')
axes[1].plot(x_data,y_data,color='blue',linewidth=2.0,linestyle='-.',label='X-Partition',marker='^')
axes[1].set_xlabel('The Size of M of Matrix', fontsize=15)
axes[1].set_ylabel('Q(I/O Complexity)', fontsize=15)
axes[1].legend(loc="best", fontsize=15)

y_data=[]
y_data2=[]
x_data=[]
max_temp=0
for i in range(38):
    m = 200 + 100 * i
    n = int(40000 / m)
    S = 40000
    x_data.append(m)

    R1S = min(m * n, S)
    R2S = min(m * n, S)

    T1S = m * m
    T2S = m * m
    T3S = m * m
    X = max(m * n, S)

    for k in range(n):
        Sre0 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
        T1S = min(math.floor(max(m * n - (k - 1) * (n - (k - 1) / 2) - S, m)), T1S)

        # Sre1 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
        Sre1 = Sre0 - (m - k + 1)
        Srem1 = Sre0 - (3 * m - 3 * k + 1)
        T21S = max(3 * (m - k) + 1 - max((math.ceil(Sre0)), 0), 0)
        T31S = max((m - k + 1) - max((math.ceil(Sre0)), 0), 0)

        Srem2 = Srem1 - (5 * m - 5 * k + 2) * (n - k + 1)
        Sre2 = Sre1 - (m - k + 1) * (n - k + 1)
        T22S = max((5 * m - 5 * k + 2) * (n - k + 1) - max(math.ceil(Srem1), 0), 0)
        T32S = max((m - k + 1) * (n - k + 1) - max(math.ceil(Sre1), 0), 0)

        Sre3 = Sre2
        Srem3 = Srem2
        T23S = max(2 * (m - k + 1) * (n - k + 1) - S, 0)
        T33S = max((m - k + 1) * (n - k + 1) - S, 0)

        T24S = max(m * n - k * (n - k / 2) - S, 0) + (m + n + 1) + max(2 * m * k - math.ceil(Srem3), 0)
        T34S = max(m * n - k * (n - k / 2) - S, 0)

        T2S = min(T21S + T22S + T23S + T24S, T2S)
        T3S = min(T31S + T32S + T33S + T34S, T3S)


    TS = 0
    if R1S - T2S <= R2S - T1S - T3S:
        TS = X - R2S + T1S + T3S
    else:
        TS = X - R1S + T2S

    TS1 = X - max(R1S, R2S) + min(T1S, T2S)

    V_1 = 4 * m * n + (4 * n - 4 * n * n) / 2
    V_2 = (3*m*n*n+5*m*n-2*n*n-n*n*n+n) / 2
    V_3 = (n - n * n * n) / 6 + (m * n + m * n * n) / 2
    V_4 = 5*n / 6 - m * n + m * n * n +n*n/2- n * n * n / 3
    V_5 = m * n * n - (n - 1) * n * (4 * n + 1) / 12
    V_all = math.ceil(V_1 + V_2 + V_3 + V_4 + V_5)

    h1 = m + 1 - max(math.floor(n + 1 - math.sqrt(max(n * n - 2 * m * n + 2 * S, 0))), 1)
    phi_1 = 4 * h1
    Phi_1 = h1

    h2 = S + Phi_1
    k2 = m + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
    i2 = n + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
    phi_2 = 3 * k2 * i2 - i2 + k2
    Phi_2 = k2 * i2

    h3 = S + Phi_1 + Phi_2
    k_3 = max(min(math.floor(
        (-math.sqrt(max(m * m + 4 * n * n + 1 + 4 * n - 8 * m * n - 4 * m + 6 * h3, 0)) + m + 2 * n + 4) / 3), n), 1)
    phi_3 = (m - k_3 + 1) * (n - k_3 + 1)
    Phi_3 = (m - k_3 + 1) * (n - k_3 + 1)

    h_4 = (k_3 - 1) * (k_3 - 1) / 2 + k_3 * (m - k_3 + 1) + (m - k_3) * (n - k_3)
    k_4 = max(min(math.floor(n - math.sqrt(max(n * n - n * m * 2 + 2 * h_4 - 1, 0))), n), 1)
    phi_4 = -k_4 * k_4 + (2 * m + 2) * k_4 - 2 * m

    V_max = S + phi_1 + phi_2 + phi_3 + phi_4
    print(m,n,S,TS/TS1)
    y_data.append(TS1 * (V_all / V_max - 1))
    y_data2.append(TS * (V_all / V_max - 1))
    max_temp=max(max_temp,TS/TS1)

print(max_temp," ")
axes[2].grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
axes[2].grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)
axes[2].set_title("(c)S=40000,M>=N,N=INT(S/M)", fontsize=15)
axes[2].plot(x_data,y_data2,color='red',linewidth=2.5,linestyle='--',label='(X1,X2)-Partition',marker='o')
axes[2].plot(x_data,y_data,color='blue',linewidth=2.0,linestyle='-.',label='X-Partition',marker='^')
axes[2].set_xlabel('The Size of M of Matrix', fontsize=15)
axes[2].set_ylabel('Q(I/O Complexity)', fontsize=15)
axes[2].legend(loc="best", fontsize=15)
plt.show()

# y_data=[]
# y_data2=[]
# x_data=[]
# max_temp=0
# for i in range(10):
#     m = 50 * (i+1)
#     n = 50 * (i+1)
#     S = m*n
#     x_data.append(m)
#
#     R1S = min(m * n, S)
#     R2S = min(m * n, S)
#
#     T1S = m * m
#     T2S = m * m
#     T3S = m * m
#     X = max(m * n, S)
#
#     for k in range(n):
#         Sre0 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
#         T1S = min(math.floor(max(m * n - (k - 1) * (n - (k - 1) / 2) - S, m)), T1S)
#
#         # Sre1 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
#         Sre1 = Sre0 - (m - k + 1)
#         Srem1 = Sre0 - (3 * m - 3 * k + 1)
#         T21S = max(3 * (m - k) + 1 - max((math.ceil(Sre0)), 0), 0)
#         T31S = max((m - k + 1) - max((math.ceil(Sre0)), 0), 0)
#
#         Srem2 = Srem1 - (5 * m - 5 * k + 1) * (n - k + 1)
#         Sre2 = Sre1 - (m - k + 1) * (n - k + 1)
#         T22S = max((5 * m - 5 * k + 1) * (n - k + 1) - max(math.ceil(Srem1), 0), 0)
#         T32S = max((m - k + 1) * (n - k + 1) - max(math.ceil(Sre1), 0), 0)
#
#         Sre3 = Sre2
#         Srem3 = Srem2
#         T23S = max(2 * (m - k + 1) * (n - k + 1) - S, 0)
#         T33S = max((m - k + 1) * (n - k + 1) - S, 0)
#
#         T24S = max(m * n - k * (n - k / 2) - Srem3, 0) + (m + n + 1) + max(2 * m * k - k - math.ceil(Srem3), 0)
#         T34S = max(m * n - k * (n - k / 2) - S, 0)
#
#         T2S = min(T21S + T22S + T23S + T24S, T2S)
#         T3S = min(T31S + T32S + T33S + T34S, T3S)
#
#     TS = 0
#     if R1S - T2S <= R2S - T1S - T3S:
#         TS = X - R2S + T1S + T3S
#     else:
#         TS = X - R1S + T2S
#
#     TS1 = X - max(R1S, R2S) + min(T1S, T2S)
#
#     V_1 = 4 * m * n + (4 * n - 4 * n * n) / 2
#     V_2 = (3*m*n*n+5*m*n-2*n*n-n*n*n+n) / 2
#     V_3 = (n - n * n * n) / 6 + (m * n + m * n * n) / 2
#     V_4 = 5*n / 6 - m * n + m * n * n +n*n/2- n * n * n / 3
#     V_5 = m * n * n - (n - 1) * n * (4 * n + 1) / 12
#     V_all = math.ceil(V_1 + V_2 + V_3 + V_4 + V_5)
#
#     h1 = m + 1 - max(math.floor(n + 1 - math.sqrt(max(n * n - 2 * m * n + 2 * S, 0))), 1)
#     phi_1 = 4 * h1
#     Phi_1 = h1
#
#     h2 = S + Phi_1
#     k2 = m + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
#     i2 = n + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
#     phi_2 = 3 * k2 * i2 - i2 + k2
#     Phi_2 = k2 * i2
#
#     h3 = S + Phi_1 + Phi_2
#     k_3 = max(min(math.floor(
#         (-math.sqrt(max(m * m + 4 * n * n + 1 + 4 * n - 8 * m * n - 4 * m + 6 * h3, 0)) + m + 2 * n + 4) / 3), n), 1)
#     phi_3 = (m - k_3 + 1) * (n - k_3 + 1)
#     Phi_3 = (m - k_3 + 1) * (n - k_3 + 1)
#
#     h_4 = (k_3 - 1) * (k_3 - 1) / 2 + k_3 * (m - k_3 + 1) + (m - k_3) * (n - k_3)
#     k_4 = max(min(math.floor(n - math.sqrt(max(n * n - n * m * 2 + 2 * h_4 - 1, 0))), n), 1)
#     phi_4 = -k_4 * k_4 + (2 * m + 2) * k_4 - 2 * m
#     V_max = S+phi_1+phi_2+phi_3+phi_4
#     print(m, n, S, TS / TS1)
#     y_data.append(TS1 * (V_all / V_max - 1))
#     y_data2.append(TS * (V_all / V_max - 1))
#     max_temp=max(TS/TS1,max_temp)
#
# print(max_temp," ")
# axes[1][0].grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
# axes[1][0].grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)
#
# axes[1][0].set_title("(d)M=N,S=MN", fontsize=15)
# axes[1][0].plot(x_data,y_data2,color='red',linewidth=2.5,linestyle='--',label='(X1,X2)-Partition',marker='o')
# axes[1][0].plot(x_data,y_data,color='blue',linewidth=2.0,linestyle='-.',label='X-Partition',marker='^')
# axes[1][0].set_xlabel('The Size of M of Matrix', fontsize=15)
# axes[1][0].set_ylabel('Q(I/O Complexity)', fontsize=15)
# axes[1][0].legend(loc="best", fontsize=15)
#
#
# y_data=[]
# y_data2=[]
# x_data=[]
# max_temp=0
# for i in range(10):
#     m = 2*50 * (i+1)
#     n = 50 * (i+1)
#     S = m*n
#     x_data.append(m)
#
#     R1S = min(m * n, S)
#     R2S = min(m * n, S)
#
#     T1S = m * m
#     T2S = m * m
#     T3S = m * m
#     X = max(m * n, S)
#
#     for k in range(n):
#         Sre0 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
#         T1S = min(math.floor(max(m * n - (k - 1) * (n - (k - 1) / 2) - S, m)), T1S)
#
#         # Sre1 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
#         Sre1 = Sre0 - (m - k + 1)
#         Srem1 = Sre0 - (3 * m - 3 * k + 1)
#         T21S = max(3 * (m - k) + 1 - max((math.ceil(Sre0)), 0), 0)
#         T31S = max((m - k + 1) - max((math.ceil(Sre0)), 0), 0)
#
#         Srem2 = Srem1 - (5 * m - 5 * k + 1) * (n - k + 1)
#         Sre2 = Sre1 - (m - k + 1) * (n - k + 1)
#         T22S = max((5 * m - 5 * k + 1) * (n - k + 1) - max(math.ceil(Srem1), 0), 0)
#         T32S = max((m - k + 1) * (n - k + 1) - max(math.ceil(Sre1), 0), 0)
#
#         Sre3 = Sre2
#         Srem3 = Srem2
#         T23S = max(2 * (m - k + 1) * (n - k + 1) - S, 0)
#         T33S = max((m - k + 1) * (n - k + 1) - S, 0)
#
#         T24S = max(m * n - k * (n - k / 2) - Srem3, 0) + (m + n + 1) + max(2 * m * k - k - math.ceil(Srem3), 0)
#         T34S = max(m * n - k * (n - k / 2) - S, 0)
#
#         T2S = min(T21S + T22S + T23S + T24S, T2S)
#         T3S = min(T31S + T32S + T33S + T34S, T3S)
#     TS = 0
#     if R1S-T2S <= R2S-T1S-T3S:
#         TS = X-R2S+T1S+T3S
#     else:
#         TS = X-R1S+T2S
#
#     TS1 = X - max(R1S,R2S) + min(T1S,T2S)
#
#     V_1 = 4 * m * n + (4 * n - 4 * n * n) / 2
#     V_2 = (3*m*n*n+5*m*n-2*n*n-n*n*n+n) / 2
#     V_3 = (n - n * n * n) / 6 + (m * n + m * n * n) / 2
#     V_4 = 5*n / 6 - m * n + m * n * n +n*n/2- n * n * n / 3
#     V_5 = m * n * n - (n - 1) * n * (4 * n + 1) / 12
#     V_all = math.ceil(V_1 + V_2 + V_3 + V_4 + V_5)
#
#     h1 = m + 1 - max(math.floor(n + 1 - math.sqrt(max(n * n - 2 * m * n + 2 * S, 0))), 1)
#     phi_1 = 4 * h1
#     Phi_1 = h1
#
#     h2 = S + Phi_1
#     k2 = m + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
#     i2 = n + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
#     phi_2 = 3 * k2 * i2 - i2 + k2
#     Phi_2 = k2 * i2
#
#     h3 = S + Phi_1 + Phi_2
#     k_3 = max(min(math.floor(
#         (-math.sqrt(max(m * m + 4 * n * n + 1 + 4 * n - 8 * m * n - 4 * m + 6 * h3, 0)) + m + 2 * n + 4) / 3), n), 1)
#     phi_3 = (m - k_3 + 1) * (n - k_3 + 1)
#     Phi_3 = (m - k_3 + 1) * (n - k_3 + 1)
#
#     h_4 = (k_3 - 1) * (k_3 - 1) / 2 + k_3 * (m - k_3 + 1) + (m - k_3) * (n - k_3)
#     k_4 = max(min(math.floor(n - math.sqrt(max(n * n - n * m * 2 + 2 * h_4 - 1, 0))), n), 1)
#     phi_4 = -k_4 * k_4 + (2 * m + 2) * k_4 - 2 * m
#
#     V_max = S + phi_1 + phi_2 + phi_3 + phi_4
#     print(m, n, S, TS / TS1)
#     y_data.append(TS1 * (V_all / V_max - 1))
#     y_data2.append(TS * (V_all / V_max - 1))
#     max_temp=max(TS/TS1,max_temp)
#
# print(max_temp," ")
# axes[1][1].grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
# axes[1][1].grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)
# axes[1][1].set_title("(e)M=2N,S=MN", fontsize=15)
# axes[1][1].plot(x_data,y_data2,color='red',linewidth=2.5,linestyle='--',label='(X1,X2)-Partition',marker='o')
# axes[1][1].plot(x_data,y_data,color='blue',linewidth=2.0,linestyle='-.',label='X-Partition',marker='^')
# axes[1][1].set_xlabel('The Size of M of Matrix', fontsize=15)
# axes[1][1].set_ylabel('Q(I/O Complexity)', fontsize=15)
# axes[1][1].legend(loc="best", fontsize=15)
#
# y_data=[]
# y_data2=[]
# x_data=[]
# max_temp=0
# for i in range(10):
#     m = 4*50 * (i+1)
#     n = 50 * (i+1)
#     S = m*n
#     x_data.append(m)
#
#     R1S = min(m * n, S)
#     R2S = min(m * n, S)
#
#     T1S = m * m
#     T2S = m * m
#     T3S = m * m
#     X = max(m * n, S)
#
#     for k in range(n):
#         Sre0 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
#         T1S = min(math.floor(max(m * n - (k - 1) * (n - (k - 1) / 2) - S, m)), T1S)
#
#         # Sre1 = S - (m * n - (k - 1) * (n - (k - 1) / 2))
#         Sre1 = Sre0 - (m - k + 1)
#         Srem1 = Sre0 - (3 * m - 3 * k + 1)
#         T21S = max(3 * (m - k) + 1 - max((math.ceil(Sre0)), 0), 0)
#         T31S = max((m - k + 1) - max((math.ceil(Sre0)), 0), 0)
#
#         Srem2 = Srem1 - (5 * m - 5 * k + 1) * (n - k + 1)
#         Sre2 = Sre1 - (m - k + 1) * (n - k + 1)
#         T22S = max((5 * m - 5 * k + 1) * (n - k + 1) - max(math.ceil(Srem1), 0), 0)
#         T32S = max((m - k + 1) * (n - k + 1) - max(math.ceil(Sre1), 0), 0)
#
#         Sre3 = Sre2
#         Srem3 = Srem2
#         T23S = max(2 * (m - k + 1) * (n - k + 1) - S, 0)
#         T33S = max((m - k + 1) * (n - k + 1) - S, 0)
#
#         T24S = max(m * n - k * (n - k / 2) - Srem3, 0) + (m + n + 1) + max(2 * m * k - k - math.ceil(Srem3), 0)
#         T34S = max(m * n - k * (n - k / 2) - S, 0)
#
#         T2S = min(T21S + T22S + T23S + T24S, T2S)
#         T3S = min(T31S + T32S + T33S + T34S, T3S)
#
#     TS = 0
#     if R1S - T2S <= R2S - T1S - T3S:
#         TS = X - R2S + T1S + T3S
#     else:
#         TS = X - R1S + T2S
#
#     TS1 = X - max(R1S, R2S) + min(T1S, T2S)
#
#     V_1 = 4 * m * n + (4 * n - 4 * n * n) / 2
#     V_2 = (3*m*n*n+5*m*n-2*n*n-n*n*n+n) / 2
#     V_3 = (n - n * n * n) / 6 + (m * n + m * n * n) / 2
#     V_4 = 5*n / 6 - m * n + m * n * n +n*n/2- n * n * n / 3
#     V_5 = m * n * n - (n - 1) * n * (4 * n + 1) / 12
#     V_all = math.ceil(V_1 + V_2 + V_3 + V_4 + V_5)
#
#     h1 = m + 1 - max(math.floor(n + 1 - math.sqrt(max(n * n - 2 * m * n + 2 * S, 0))), 1)
#     phi_1 = 4 * h1
#     Phi_1 = h1
#
#     h2 = S + Phi_1
#     k2 = m + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
#     i2 = n + 1 - max(math.floor(n + 2 - math.sqrt(max((n + 1) * (n + 1) - 2 * m * n + 2 * h1 - 2 * m, 0))), 1)
#     phi_2 = 3 * k2 * i2 - i2 + k2
#     Phi_2 = k2 * i2
#
#     h3 = S + Phi_1 + Phi_2
#     k_3 = max(min(math.floor(
#         (-math.sqrt(max(m * m + 4 * n * n + 1 + 4 * n - 8 * m * n - 4 * m + 6 * h3, 0)) + m + 2 * n + 4) / 3), n), 1)
#     phi_3 = (m - k_3 + 1) * (n - k_3 + 1)
#     Phi_3 = (m - k_3 + 1) * (n - k_3 + 1)
#
#     h_4 = (k_3 - 1) * (k_3 - 1) / 2 + k_3 * (m - k_3 + 1) + (m - k_3) * (n - k_3)
#     k_4 = max(min(math.floor(n - math.sqrt(max(n * n - n * m * 2 + 2 * h_4 - 1, 0))), n), 1)
#     phi_4 = -k_4 * k_4 + (2 * m + 2) * k_4 - 2 * m
#
#     V_max = S + phi_1 + phi_2 + phi_3 + phi_4
#     print(m, n, S, TS / TS1)
#     y_data.append(TS1 * (V_all / V_max - 1))
#     y_data2.append(TS * (V_all / V_max - 1))
#     max_temp=max(TS/TS1,max_temp)
#
# print(max_temp," ")
# print("111")
# axes[1][2].grid(which="major", axis='x', color='#DAD8D7', alpha=0.5, zorder=1)
# axes[1][2].grid(which="major", axis='y', color='#DAD8D7', alpha=0.5, zorder=1)
#
# axes[1][2].plot(x_data,y_data2,color='red',linewidth=2.5,linestyle='--',label='(X1,X2)-Partition',marker='o')
# axes[1][2].plot(x_data,y_data,color='blue',linewidth=2.0,linestyle='-.',label='X-Partition',marker='^')
# axes[1][2].set_xlabel('The Size of M of Matrix', fontsize=15)
# axes[1][2].set_ylabel('Q(I/O Complexity)', fontsize=15)
# axes[1][2].set_title("(f)M=4N,S=MN", fontsize=15)
# axes[1][2].legend(loc="best", fontsize=15)