import sympy as sp

# h-tiered clos topology
# number of variables = h*(h+1)/2
# number of equations = 2*h
# X[i][j] is the variable which need to be solved.
# t[i][j] is not a variable in the equation, it is a fixed value.
# C is a constant value.

# 2*X[1][1] + 2*X[1][2] + ... + 2*X[1][h-1] + X[1][h] = 1
# .......
# 2*X[h-2][1] + 2*X[h-2][2] + X[h-2][3] = 1
# 2*X[h-1][1] + X[h-1][2] = 1
# X[h][1] = 1

# t[1][1]*X[1][h] = C
# t[2][1]*X[1][h-1] + t[2][2]*X[2][h-1] = C
# t[3][1]*X[1][h-2] + t[3][2]*X[2][h-2] + t[3][3]*X[3][h-2] = C
# ......
# t[h][1]*X[1][1] + t[h][2]*X[2][1] + ... + t[h][h]*X[h][1] = C


# h = 2
# t11, t21, t22 = sp.symbols('t11 t21 t22')
# t11, t21, t22 = 28, 28, 1
# C = sp.symbols('C')
# A = sp.Matrix([
#     [2, 1, 0],
#     [0, 0, 1],
#     [0, t11, 0],
#     [t21, 0, t22]
# ])
# B = sp.Matrix(
#     [1, 1, C, C]
# )

# h = 3 
# t11, t21, t22, t31, t32, t33 = sp.symbols('t11 t21 t22 t31 t32 t33')
# C = sp.symbols('C')
# A = sp.Matrix([
#     [2, 2, 1, 0, 0, 0],
#     [0, 0, 0, 2, 1, 0],
#     [0, 0, 0, 0, 0, 1],
#     [0, 0, t11, 0, 0, 0],
#     [0, t21, 0, 0, t22, 0],
#     [t31, 0, 0, t32, 0, t33]
# ])
# B = sp.Matrix(
#     [1, 1, 1, C, C, C]
# )

# h = 4
# t11, t21, t22, t31, t32, t33, t41, t42, t43, t44 = sp.symbols('t11 t21 t22 t31 t32 t33 t41 t42 t43 t44')
# C = sp.symbols('C')
# A = sp.Matrix([
#     [2, 2, 2, 1, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 2, 2, 1, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 2, 1, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
#     [0, 0, 0, t11, 0, 0, 0, 0, 0, 0],
#     [0, 0, t21, 0, 0, 0, t22, 0, 0, 0],
#     [0, t31, 0, 0, 0, t32, 0, 0, t33, 0],
#     [t41, 0, 0, 0, t42, 0, 0, t43, 0, t44]
# ])
# B = sp.Matrix(
#     [1, 1, 1, 1, C, C, C, C]
# )


# if A.shape[0] == A.shape[1]:
#     # 计算系数矩阵 A 的行列式
#     det_A = A.det()
#     # 打印行列式
#     print(f"系数矩阵 A 的行列式为: {det_A}")

# # 计算系数矩阵 A 的秩和增广矩阵 [A|B] 的秩
# # print(A, B)
# rank_A = A.rank()
# augmented_matrix = A.row_join(B)
# rank_augmented = augmented_matrix.rank()

# # 打印秩
# print(f"系数矩阵 A 的秩为: {rank_A}")
# print(f"增广矩阵 [A|B] 的秩为: {rank_augmented}")

# # 判断解的存在性
# if rank_A == rank_augmented:
#     if rank_A == A.shape[1]:
#         print("秩等于变量数，方程组有唯一解。")
#     else:
#         print("秩等于增广矩阵的秩，但小于变量数，方程组有无穷多解。")
# else:
#     print("系数矩阵的秩小于增广矩阵的秩，方程组无解。")

def generate_matrices(h):
    # 符号变量生成
    num = h*(h+1)//2
    t = sp.symbols(f't:{num}')
    C = sp.symbols('C')
    
    # 生成系数矩阵 A
    A = []
    tp = 0
    for i in range(h, 1, -1):
        row = [0] * num
        for j in range(tp, tp+i):
            row[j] = 2
        row[tp+i-1] = 1
        tp += i
        A.append(row)

    tp, idx = h-1, 0
    for i in range(1, h+1):
        row = [0] * num
        st = tp
        for j in range(i):
            row[st] = t[idx]
            idx += 1
            st += h-j
        tp -= 1
        A.append(row)
    
    A = sp.Matrix(A)
    print(A)
    
    # 生成常数向量 B
    B = sp.Matrix([1] * (h-1) + [C] * h)
    print(B)
    
    return A, B

def is_solvable(h):
    A, B = generate_matrices(h)
    
    # 计算系数矩阵 A 的秩
    rank_A = A.rank()
    # 生成增广矩阵 [A|B]
    augmented_matrix = A.row_join(B)
    # 计算增广矩阵 [A|B] 的秩
    rank_augmented = augmented_matrix.rank()

    # 使用高斯消元法化简矩阵
    rref_matrix, pivot_columns = A.rref()
    print("系数矩阵 A 化简后的矩阵 (RREF):")
    sp.pprint(rref_matrix)
    rref_matrix, pivot_columns = augmented_matrix.rref()
    print("增广矩阵 [A|B] 化简后的矩阵 (RREF):")
    sp.pprint(rref_matrix)

    print(f"系数矩阵 A 的秩为: {rank_A}, 增广矩阵 [A|B] 的秩为: {rank_augmented}")
    
    # 判断解的存在性
    if rank_A == rank_augmented:
        if rank_A == A.shape[1]:
            print("秩等于变量数，方程组有唯一解。")
        else:
            print("秩等于增广矩阵的秩，但小于变量数，方程组有无穷多解。")
    else:
        print("系数矩阵的秩小于增广矩阵的秩，方程组无解。")

if __name__ == "__main__":
    h = 2
    for h in range(2, 100):
        print(f"==================== h = {h} ====================")
        is_solvable(h)

# 1-hop 1
# 3-hop x y x (2*x + y = 1)
# 5-hop a b c a b (2*a + 2*b + c = 1)
# ====== 5介数中心性 * c
# Core: (k**4-k**3)//8 * c
# ====== 5介数中心性 * b + 3介数中心性 * y
# Agg: (k**4-k**3)//8 * b + (k**2-2*k)//8 * y
# ====== 5介数中心性 * a + 3介数中心性 * x + 1
# Edge: (k**4-k**3)//8 * a + (k**2-2*k)//4 * x + 1
# c = (16+k**4-k**3+2*k**2-4*k) / 5*(k**4-k**3)
# a = (k**4-k**3-k**2+k-4) / 5*(k**4-k**3)
# b = (k**4-k**3+k-4) / 5*(k**4-k**3)
# x = (3*k**2-5*k-20) / 10*k*(k-2)
# y = (2*k**2-5*k+20) / 5*k*(k-2)

# spine leaf
# 1-hop 1
# 3-hop x y x (2*x + y = 1)
# ====== 3介数中心性 * y
# Spine: m*(m-1)//2 * y 
# ====== 3介数中心性 * x + 1
# Leaf: n*(m-1) * x + 1
# x = (m*(m-1)-2) // 2*(m-1)*(m+n)
# y = (2*(m-1)*n+4) // 2*(m-1)*(m+n) 