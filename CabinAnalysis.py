import numpy as np
import Parameters as pm
import trusspy as tp
import matplotlib.pyplot as plt


concept = pm.ConceptParameters(0)
h = concept.cabin.H_cabin
l = concept.cabin.L_cabin
w = concept.cabin.W_cabin
M = tp.Model()

N1 = tp.Node(1, (0,0,0))
N2 = tp.Node(2,(0,0,h))
N3 = tp.Node(3,(0,w,h))
N4 = tp.Node(4,(0,w,0))
N5 = tp.Node(5,(0,0.75*w,0))
N6 = tp.Node(6,(0,0.25*w,0))
N7 = tp.Node(7,(l,0,h))
N8 = tp.Node(8,(l,w,h))
N9 = tp.Node(9,(l,w,0))
N10 = tp.Node(10,(l,0,0))

E_mod = 1
A_truss = 1

EA = tp.Element(1,[1,2],material_properties=[E_mod],geometric_properties=[A_truss])
EB = tp.Element(2,[2,6],material_properties=[E_mod],geometric_properties=[A_truss])
EC = tp.Element(3,[2,3],material_properties=[E_mod],geometric_properties=[A_truss])
ED = tp.Element(4,[3,5],material_properties=[E_mod],geometric_properties=[A_truss])
EE = tp.Element(5,[3,4],material_properties=[E_mod],geometric_properties=[A_truss])
EF = tp.Element(6,[4,5],material_properties=[E_mod],geometric_properties=[A_truss])
EG = tp.Element(7,[5,6],material_properties=[E_mod],geometric_properties=[A_truss])
EH = tp.Element(8,[6,1],material_properties=[E_mod],geometric_properties=[A_truss])
EI = tp.Element(9,[2,7],material_properties=[E_mod],geometric_properties=[A_truss])
EJ = tp.Element(10,[7,8],material_properties=[E_mod],geometric_properties=[A_truss])
EK = tp.Element(11,[8,3],material_properties=[E_mod],geometric_properties=[A_truss])
EL = tp.Element(12,[8,9],material_properties=[E_mod],geometric_properties=[A_truss])
EM = tp.Element(13,[9,10],material_properties=[E_mod],geometric_properties=[A_truss])
EN = tp.Element(14,[10,7],material_properties=[E_mod],geometric_properties=[A_truss])
EO = tp.Element(15,[9,4],material_properties=[E_mod],geometric_properties=[A_truss])
EP = tp.Element(16,[9,5],material_properties=[E_mod],geometric_properties=[A_truss])
EQ = tp.Element(17,[6,10],material_properties=[E_mod],geometric_properties=[A_truss])
ER = tp.Element(18,[1,10],material_properties=[E_mod],geometric_properties=[A_truss])


M.Nodes.add_nodes([N1,N2,N3,N4,N5,N6,N7,N8,N9,N10])
M.Elements.add_elements([EA,EB,EC,ED,EE,EF,EG,EH,EI,EJ,EK,EL,EM,EN,EO,EP,EQ,ER])



with M.Boundaries as MB:
    MB.add_bound_U(5, (0,0,0))
    MB.add_bound_U(6, (0,0,0))
    MB.add_bound_U(9, (0,0,1))
    MB.add_bound_U(4, (0,0,1))
    MB.add_bound_U(2, (0,0,1))
    MB.add_bound_U(7, (0,0,1))
    MB.add_bound_U(8, (0,0,1))
    MB.add_bound_U(3, (0,0,1))

with M.ExtForces as ME:
    ME.add_force(9,(0,0,1))
    ME.add_force(4,(0,0,1))
    ME.add_force(2,(0,0,1))
    ME.add_force(7,(0,0,1))
    ME.add_force(8,(0,0,1))
    ME.add_force(3,(0,0,1))


# F1 = tp.ExternalForce(1, (1,1,1))
# F2 = tp.ExternalForce(10, (1,1,1))
# F3 = tp.ExternalForce(9, (1,1,1))
# F4 = tp.ExternalForce(4, (1,1,1))



# M.Boundaries.add_bounds_U([B1,B2,B3,B4])
# M.ExtForces.add_forces([F1,F2,F3,F4])

M.build()
M.run()

fig, ax = M.plot_model(config=['deformed'],
                       view='3d',
                       contour='force',
                       force_scale=500.0,
                       inc=40)
plt.show()