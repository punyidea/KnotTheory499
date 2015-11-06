from smith_normal import smith_normal_form
import numpy as np

mat_sizes =  np.random.randint(3,50,(40,2))

for mat_size in mat_sizes:


    #A = np.array([[8,0,0,0,0],[0,4,6,18,-1],[0,7,6,21,0],[0,0,0,0,0],[0,3,2,3,1]])
    A = np.random.randint(-20,20,mat_size)
    B = smith_normal_form(A)

    print(A)
    print(B)

