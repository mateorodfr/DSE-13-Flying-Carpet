import numpy as np
import Parameters as pm






#Defining box and node points (where the rod bois connect)
B_width = 0.5
W , H , D = 1, 1, 1
Point_1 = np.array([D/2,W/2,0])
Point_2 = np.array([-D/2,W/2,0])
Point_3 = np.array([-D/2,-W/2,0])
Point_4 = np.array([D/2,-W/2,0])
Point_5 = np.array([D/2,W/2,H])
Point_6 = np.array([-D/2,W/2,H])
Point_7 = np.array([-D/2,-W/2,H])
Point_8 = np.array([D/2,-W/2,H])
Point_B1 = np.array([D/2,B_width/2,H])
Point_B2 = np.array([D/2,-B_width/2,H])

Pointlist = np.array([Point_1,Point_2,Point_3,Point_4,Point_5,Point_6,Point_7,Point_8,Point_B1,Point_B2])

#Defining boundary conditions (forces at the node bois)

F_1 = np.array([1,1,1])
F_2 = np.array([1,1,1])
F_3 = np.array([1,1,1])
F_4 = np.array([1,1,1])
F_5 = np.array([1,1,1])
F_6 = np.array([1,1,1])
F_7 = np.array([1,1,1])
F_8 = np.array([1,1,1])
F_B1 = np.array([1,1,1])
F_B2 = np.array([1,1,1])

Forcelist = np.array([F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_B1,F_B2])

#Defining connector function

def connector(Point_1,Point_2):  # IMPORTANT rod goes from point 1 to point 2 (2 is positive)
    vec = Point_2-Point_1
    length = np.sqrt( (Point_2[0]-Point_1[0])**2+(Point_2[1]-Point_1[1])**2+(Point_2[2]-Point_1[2])**2)
    unitvec = vec / length 
    Rod = np.array([Point_1 , Point_2 , vec , length , unitvec])
    return 

#Making the rods we want to test, below lies a perfect box of sides

Rod1 = connector(Point_1,Point_2)
Rod2 = connector(Point_2,Point_3)
Rod3 = connector(Point_3,Point_4)
Rod4 = connector(Point_4,Point_1)
Rod5 = connector(Point_5,Point_6)
Rod6 = connector(Point_6,Point_7)
Rod7 = connector(Point_7,Point_8)
Rod8 = connector(Point_8,Point_B2)
Rod9 = connector(Point_B2,Point_B1)
Rod10 = connector(Point_B1,Point_5)
Rod11 = connector(Point_5,Point_1)
Rod12 = connector(Point_8,Point_4)
Rod13 = connector(Point_7,Point_3)
Rod14 = connector(Point_6,Point_2)

#Here I add those rods which are necessary but don't lead to a perfect box

Rod15 = connector(Point_B1,Point_1)
Rod16 = connector(Point_B2,Point_4)
Rod17 = connector(Point_8,Point_3)
Rod18 = connector(Point_5,Point_2)
Rod19 = connector(Point_6,Point_3)

#Making Rod array for computing

Rodarray = np.array([Rod1,Rod2,Rod3,Rod4,Rod5,Rod6,Rod7,Rod8,Rod9,Rod10,Rod11,Rod12,Rod13,Rod14,Rod15,Rod16,Rod17,Rod18,Rod19])



#Some sorting tests
print(Rodarray)




####Matrix operations gonna be below



#Nodenum = 10
#Matrixsize = np.zeros( 3 * Nodenum , 3 * Nodenum )
