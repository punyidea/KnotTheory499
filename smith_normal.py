__author__ = 'vsp'
import numpy as np
def smith_normal_form(A):
    S = A.copy()
    nrows, ncols = np.shape(S)
    piv_on_row =  nrows <= ncols

    if not piv_on_row:
        nrows,ncols = ncols,nrows
        S = S.T

    def reduce_rows(piv,other_piv):
        nonlocal S
        #reduce all rows
        for other_row_ind,other_val in enumerate(S[piv+1:, other_piv]):
            if other_val:
                other_row_ind = other_row_ind+piv+1
                piv_val = S[piv, other_piv]
                piv_row = S[piv, :]

                if other_val % piv_val:
                    common_divisor,sigma,tau = extend_euclid(piv_val,other_val)
                    alpha = piv_val//common_divisor
                    beta = other_val//common_divisor

                    Lo = np.eye(nrows,dtype=int)
                    Lo[[piv,piv,other_row_ind,other_row_ind],
                       [piv,other_row_ind,piv,other_row_ind]] = \
                        [sigma,tau,-beta,alpha]

                    S = np.dot(Lo, S)

                else:
                    #reduce the row by the operation
                    S[other_row_ind, :] -= other_val // piv_val * piv_row

        return

    def reduce_cols(piv,other_piv):
        nonlocal S

        #reduce all columns, given
        for other_col_ind,other_val in enumerate(S[piv, other_piv+1:]):
            if other_val:
                other_col_ind = other_col_ind+other_piv+1
                piv_val = S[piv, other_piv]
                piv_col = S[:, other_piv]

                if other_val % piv_val:
                    common_divisor,sigma,tau = extend_euclid(piv_val,other_val)
                    alpha = piv_val//common_divisor
                    beta = other_val//common_divisor

                    Lo = np.eye(ncols,dtype=int)
                    Lo[[piv,piv,other_col_ind,other_col_ind],
                       [piv,other_col_ind,piv,other_col_ind]] = \
                        [sigma,-beta,tau,alpha]

                    S = np.dot(S, Lo)

                else:
                    #reduce the column by the operation
                    S[:, other_col_ind] -= other_val // piv_val * piv_col
        return(S)

    def reduce_piv(piv,other_piv):
        nonlocal S
        reduce_rows(piv,other_piv)
        reduce_cols(piv,other_piv)


        while any(S[piv+1:, other_piv]) or \
            any(S[piv, other_piv+1:]):

            reduce_rows(piv,other_piv)
            if any(S[piv, other_piv+1:]):
                reduce_cols(piv,other_piv)


    def step_1():
        '''
        Performs the first step of the SNF algorithm-
            diagonalization.
        :param other_piv:
        :return:
        '''
        other_piv = 0
        nonlocal S
        for piv in range(nrows):
            ideal_pivot = min_nonzero_index(S[:, other_piv])
            while ideal_pivot is None:
                other_piv+=1
                if other_piv == ncols: # Then only columns of zeroes remain.
                    return
                ideal_pivot = min_nonzero_index(S[:, other_piv])

            S[[ideal_pivot, piv], :] = S[[piv, ideal_pivot], :] #swap rows to make the best pivot

            reduce_piv(piv,other_piv)

            other_piv +=1
            if other_piv == ncols: # Then only columns of zeroes remain.
                return

        return
    def fix_cols():
        '''
        Orders the columns such that all zero-columns are the right-most columns.
        :return:
        '''
        nonlocal S
        has_all_zeros = np.logical_not(np.any(S,axis = 0)) #find which columns have only zeros
        col_indices = np.argsort(has_all_zeros,kind = 'mergesort')
        S = S[:,col_indices]


    def fix_divisibility():
        '''
        Once diagonalizing the matrix, perform operations to ensure the
        divisibility conditions of SNF.
        :return:
        '''
        nonlocal S
        for piv in range(nrows-1):
            if S[piv,piv]:
                if S[piv+1,piv+1] % S[piv,piv]:
                    S[:,piv] += S[:,piv+1]
                    reduce_piv(piv,piv)
            else: return


    step_1()
    fix_cols()
    fix_divisibility()

    if not piv_on_row:
        S = S.T

    return S





def vect_gcd(vect):
    '''
    computes the gcd of all numbers in the vector
    :param vect:
    :return:
    '''
    def vect_scal_gcd(scal,vect):
        if len(vect)==1:
            return gcd(scal,vect[0])
        else:
            return gcd(scal,vect_scal_gcd(vect[0],vect[1:]))
    if len(vect)==1:
        return vect[0]
    else:
        return vect_scal_gcd(vect[0],vect[1:])



def gcd(a,b):
    '''
    returns the greatest common denominator of two numbers.
    :param a:
    :param b:
    :return:
    '''
    a,b = int(a), int(b)
    while b !=0:
        a,b = b, a%b
    return a

def extend_euclid(a,b):
    '''
    returns the gcd of two numbers,a and b, and
    two integers ,s, and t, which have the property that
    a*s + b*t = gcd(a,b)
    :param a:
    :param b:
    :return: gcd, s, t
    '''
    a,b = int(a), int(b)
    s_past, s_now =1,0
    t_past,t_now = 0,1

    while b !=0:
        q = a//b
        a,b = b, a%b
        s_past,s_now = s_now,s_past-q*s_now
        t_past,t_now = t_now,t_past-q*t_now
    gcd,s,t = a,s_past,t_past
    return gcd,s,t

def min_nonzero_index(vect):
    min_index= min_val = None
    for ind,i in enumerate(vect):
        if i !=0:
            if not min_val:
                min_val,min_index = abs(i),ind
            elif abs(i)<min_val:
                min_val,min_index = abs(i),ind
    return min_index

