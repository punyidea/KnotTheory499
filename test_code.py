from smith_normal import smith_normal_form
import numpy as np



A = np.array([[8,0,0,0,0],[0,4,6,18,-1],[0,7,6,21,0],[0,0,0,0,0],[0,3,2,3,1]])

B = smith_normal_form(A)
C =3