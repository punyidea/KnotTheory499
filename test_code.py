from smith_normal import smith_normal_form
import numpy as np

def matrix_mathematica_print(mat):

    mat_string = '}, {'.join([', '.join(mat_row) for mat_row in mat])

    print('\{\{{}\}\}'.format(mat_string))

mat_sizes =  np.random.randint(3,45,(400,1))
mat_sizes = np.concatenate((mat_sizes,mat_sizes),axis=1)

# A = np.array([[  0,   3,  -5,  -6,  10, -16, -11,   0,  -6, -14,  -2, -17,   8,  -5,   6,   6],
#               [ 11,  -7,   9, -19,   8,  -7,   1,  15,   6,  13,   3,   7,  -5,  13,  11,  10],
#               [ 16,  -9, -10,  -7, -19,  10,  -9,   6,  -1, -11,  -9,   1,   8, -12,  14,  11],
#               [ 12,   4,  -7,   2,  18,  -6,   1,  10,   1,  12,   8,   4,  -9,  10,  13,  -8],
#               [  4,  12,  -7,  -6,  -7,  15,   6, -13,  -1, -16,  13,  10, -17,  -5,   5,   8],
#               [ -7, -16,  -1,  -8,   1,   1,   0,   6,  13,  -5,  -4,  -5,  -9, -18,  12,   5],
#               [ -2,   6, -17,  -7,  -1, -17,   8, -13,  -2,  13, -10,  15, -13,  -9,   2, -10],
#               [-13, -16,   6,  15, -20,  -6,   9, -16,  15,  12, -19,  -1,   1, -20, -13,  -2],
#               [ 15,   9, -14,  13,  -3,  11, -14,   7,   0,   0, -19,   9,  -8,  -1,   8,   9],
#               [  8, -20, -17,  -8,   8, -13,  10, -18, -20,   3,  15,  10, -11, -14, -10, -17],
#               [ 15,  18,  -5,  19,  -2,   5,  -5,  -3,   3,  -6, -12,  12,   4,  11,  -3,   9],
#               [ 11,  15,  -3,  -8, -12,  14,   6,  19,  13,  -9,   9, -15, -19,  16, -18, -14],
#               [  0,   9,   5,  13,   5, -20, -20,  -7,   5,   0,  -5, -13,  -8,  -6,   1,  -1],
#               [ -5,   0,   1, -14,   3,  -5,  15, -15,  15, -19,   0,  -4,  15,   6,  -2,   8],
#               [  4,  19, -19,  -9,   1,  -5,  10,   2,   1,   4,   1,   5,   1,   3,  10,  12],
#               [-18,   6,   5,  12,   4,  11, -15, -14,  -8, -20,   8,  10,   2, -16,   9,  -8],
#               [  6,  -8,  -4,  14, -12,  -9,  11,  -3, -19,   5, -20,  13,  -7, -17, -10, -11]]
# )
#
# B = smith_normal_form(A)

for mat_size in mat_sizes:


    #A = np.array([[8,0,0,0,0],[0,4,6,18,-1],[0,7,6,21,0],[0,0,0,0,0],[0,3,2,3,1]])
    A = np.random.randint(-2,2,mat_size)
    A = A.astype(int,casting='unsafe')
    B = smith_normal_form(A)

    if abs(round(np.linalg.det(A))) != abs(round(np.linalg.det(B))):
        #something has gone wrong in computation.
        matrix_mathematica_print(A)
        print(A)
        print(B)
        print(round(np.linalg.det(A)))
        c = 1

