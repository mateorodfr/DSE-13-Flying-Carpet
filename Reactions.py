from math import sin, cos, pi, sqrt
from scipy.linalg import solve
import numpy as np
import matplotlib.pyplot as plt


class Centerpiece:
    def __init__(self, l_a1, l_a2, l_a3, l_a4, cptype):
        """
        This method is automatically called when an instance of the Reactions object is created
        i.e. var = Reactions() creates an instance

        The 'self' object allows variables and methods to be accessed and modified from outside the class, i.e.:
        Accessing a method:     var.set_thrust(100, 100, 100, 100)
        Modifying a variable:   var.p1 = 200

        float l_a1 : length of arm 1 from hinge to point of force action
        float l_a2 : length of arm 2 from hinge to point of force action
        float l_a3 : length of arm 3 from hinge to point of force action
        float l_a4 : length of arm 4 from hinge to point of force action
        string cptype : type of centerpiece structural solution
            "beams_square" : 4 beams in a square pattern
        """

        # ----------------------
        # Initializing constants
        # ----------------------
        self.__l_a1 = l_a1
        self.__l_a2 = l_a2
        self.__l_a3 = l_a3
        self.__l_a4 = l_a4
        self.__cptype = cptype

        # ----------------------
        # Initializing variables
        # ----------------------

        # They are empty but may be used later
        # If they are used in calculations before they are set, an error is returned
        # This is nice since it prevents the thrusts from being used before they are specified explicitly

        # The double underscore before each variable name is a python hack
        # It emulates private variables - variables that should not be accessed or modified externally
        # It prevents them from being changed if you don't explicitly call the method to change them
        # i.e. the thrusts can only be changed by using name.set_p(), trying to do name.__p1 = 1 will give an error

        # Initializing engine thrust variables
        self.__p1 = None  # p1 is (from entrance perspective) the front left engine
        self.__p2 = None  # p2 is (from entrance perspective) the back left engine
        self.__p3 = None  # p3 is (from entrance perspective) the back right engine
        self.__p4 = None  # p4 is (from entrance perspective) the front right engine

        # Initializing engine torque variables
        self.__t1 = None  # t1 is (from entrance perspective) the front left engine
        self.__t2 = None  # t2 is (from entrance perspective) the front left engine
        self.__t3 = None  # t3 is (from entrance perspective) the front left engine
        self.__t4 = None  # t4 is (from entrance perspective) the front left engine

        # Initializing arm angles
        self.__th1 = None  # th1 is (from entrance perspective) the front left arm
        self.__th2 = None  # th2 is (from entrance perspective) the front left arm
        self.__th3 = None  # th3 is (from entrance perspective) the front left arm
        self.__th4 = None  # th4 is (from entrance perspective) the front left arm

        # Initialize beam properties if beams_square cptype is used
        if cptype == "beams_square":
            # Width and length of centerpiece
            self.__l1 = None  # Width
            self.__l2 = None  # Length

            self.__h = None
            self.__b = None
            self.__tf = None
            self.__tw = None

            self.__sigma_yield = None
            self.__mass = None

            # Area Moments of Inertia of beams
            self.__Ixx = None
            self.__Iyy = None

            # Wrote the program to have variable MOIs for each beam but in the end they all use the same one
            self.__Ixx_12 = None
            self.__Ixx_23 = None
            self.__Ixx_34 = None
            self.__Ixx_41 = None
            self.__Iyy_12 = None
            self.__Iyy_23 = None
            self.__Iyy_34 = None
            self.__Iyy_41 = None

            # Young's Modulus
            self.__e_modulus = None

            # Reactions at connectors
            self.__F_c1_x = None
            self.__F_c1_y = None
            self.__F_c1_z = None

            self.__F_c2_x = None
            self.__F_c2_y = None
            self.__F_c2_z = None

            self.__F_c3_x = None
            self.__F_c3_y = None
            self.__F_c3_z = None

            self.__F_c4_x = None
            self.__F_c4_y = None
            self.__F_c4_z = None

            self.__M_c1_x = None
            self.__M_c1_y = None
            self.__M_c1_z = None

            self.__M_c2_x = None
            self.__M_c2_y = None
            self.__M_c2_z = None

            self.__M_c3_x = None
            self.__M_c3_y = None
            self.__M_c3_z = None

            self.__M_c4_x = None
            self.__M_c4_y = None
            self.__M_c4_z = None

            # Internal moments
            self.__Mz_12 = None
            self.__Mz_23 = None
            self.__Mz_34 = None
            self.__Mz_41 = None

    def set_p(self, p1, p2, p3, p4):
        """
        Sets the (net) force produced at each arm tip
        This includes the weight of the engine
        Downward positive (positive z-direction)
        """
        self.__p1 = p1
        self.__p2 = p2
        self.__p3 = p3
        self.__p4 = p4

    def set_t(self, t1, t2, t3, t4):
        """
        Sets the (net) torque at each arm tip
        Positive according to right hand rule in z-direction
        """
        self.__t1 = t1
        self.__t2 = t2
        self.__t3 = t3
        self.__t4 = t4

    def set_arm_angle(self, th1, th2, th3, th4):
        """
        Sets the angles of the engine arms
        Angles are relative to x-axis
        Positive according to right hand rule along z-axis
        """
        self.__th1 = th1
        self.__th2 = th2
        self.__th3 = th3
        self.__th4 = th4

    def set_beam_config(self, l1, l2, h, b, tw, tf, e_modulus, rho, sigma_yield):
        """
        float l1 : y-distance between beams 12 and 34
        float l2 : x-distance between beams 23 and 41
        float h : web height
        float b : flange width
        float tw : web thickness
        float tf : flange thickness
        float e_modulus : E-modulus of material (same material for all beams)
        """
        if self.__cptype == "beams_square":
            self.__l1 = l1
            self.__l2 = l2

            self.__h = h
            self.__b = b
            self.__tw = tw
            self.__tf = tf

            self.__sigma_yield = sigma_yield

            A = 2 * b * tf + h * tw
            self.__mass = (A * l1 * 2 + A * l2 * 2) * rho

            # print("material: Al 6061-T6")
            # print("centerpiece mass: " + str(self.__mass) + " kg")

            # Calculate MOIs
            Ixx = h * h * h * tw / 12 + 2 * (tf * tf * tf * b / 12 + tf * b * (h + tf) * (h + tf) / 4)
            Iyy = tw * tw * tw * h / 12 + 2 * (b * b * b * tf / 12)

            # Check for validity of no torsion assumption (bending stiffness >> torsional stiffness)
            J = 1 / 3 * (2 * b * tf * tf * tf + h * tw * tw * tw)
            if round(J / Ixx * 100, 5) > 1:
                print("WARNING: Torsional stiffness significant relative to bending stiffness. "
                      "Results may be incorrect")

            # Wrote the program to have variable MOIs for each beam but in the end they all use the same one
            self.__Ixx = Ixx
            self.__Iyy = Iyy
            self.__Ixx_12 = Ixx
            self.__Ixx_23 = Ixx
            self.__Ixx_34 = Ixx
            self.__Ixx_41 = Ixx
            self.__Iyy_12 = Iyy
            self.__Iyy_23 = Iyy
            self.__Iyy_34 = Iyy
            self.__Iyy_41 = Iyy
            self.__e_modulus = e_modulus
        else:
            print("Error: set_beam_config called while cptype not a beam type")
            exit()

    def solve_c(self, ver=False):
        """
        Solves for reactions (and structure rotations and internal moments) at connectors
        """

        # Assume thrust always acts along z-axis : rotation of motor assumed negligible
        # Assume bending stiffness >> torsional stiffness
        self.__F_c1_z = self.__p1
        self.__M_c1_x = self.__p1 * self.__l_a1 * sin(self.__th1)
        self.__M_c1_y = - self.__p1 * self.__l_a1 * cos(self.__th1)
        self.__M_c1_z = self.__t1

        self.__F_c2_z = self.__p2
        self.__M_c2_x = self.__p2 * self.__l_a2 * sin(self.__th2)
        self.__M_c2_y = - self.__p2 * self.__l_a2 * cos(self.__th2)
        self.__M_c2_z = self.__t2

        self.__F_c3_z = self.__p3
        self.__M_c3_x = self.__p3 * self.__l_a3 * sin(self.__th3)
        self.__M_c3_y = - self.__p3 * self.__l_a3 * cos(self.__th3)
        self.__M_c3_z = self.__t3

        self.__F_c4_z = self.__p4
        self.__M_c4_x = self.__p4 * self.__l_a4 * sin(self.__th4)
        self.__M_c4_y = - self.__p4 * self.__l_a4 * cos(self.__th4)
        self.__M_c4_z = self.__t4

        # Check that forces were given
        if self.__t1 is None:
            print("Error: Forces not specified")
            exit()

        if self.__cptype == "beams_square":
            # Assemble matrix
            # Calculate some terms beforehand to keep matrix clean
            L1 = self.__l1
            L2 = self.__l2

            A9_3 = L2 * L2 * L2 / (3 * self.__e_modulus * self.__Iyy_12)
            A9_10 = L2 * L2 / (2 * self.__e_modulus * self.__Iyy_12)
            A10_6 = -L1 * L1 * L1 / (3 * self.__e_modulus * self.__Iyy_23)
            A10_11 = L1 * L1 / (2 * self.__e_modulus * self.__Iyy_23)
            A11_7 = - L2 * L2 * L2 / (3 * self.__e_modulus * self.__Iyy_34)
            A11_12 = L2 * L2 / (2 * self.__e_modulus * self.__Iyy_34)
            A12_2 = L1 * L1 * L1 / (3 * self.__e_modulus * self.__Iyy_41)
            A12_9 = L1 * L1 / (2 * self.__e_modulus * self.__Iyy_41)
            A13_3 = L2 * L2 / (2 * self.__e_modulus * self.__Iyy_12)
            A13_10 = L2 / (self.__e_modulus * self.__Iyy_12)
            A14_6 = - L1 * L1 / (2 * self.__e_modulus * self.__Iyy_23)
            A14_11 = L1 / (self.__e_modulus * self.__Iyy_23)
            A15_7 = - L2 * L2 / (2 * self.__e_modulus * self.__Iyy_34)
            A15_12 = L2 / (self.__e_modulus * self.__Iyy_34)
            A16_2 = L1 * L1 / (2 * self.__e_modulus * self.__Iyy_41)
            A16_9 = L1 / (self.__e_modulus * self.__Iyy_41)

            # Matrix assembly
            A = [
                [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # eq1
                [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # eq2
                [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # eq3
                [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # eq4
                [0, 0, L2, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0],  # eq5
                [0, 0, 0, 0, 0, -L1, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0],  # eq6
                [0, 0, 0, 0, 0, 0, -L2, 0, 0, 0, -1, 1, 0, 0, 0, 0],  # eq7
                [0, L1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0],  # eq8
                [0, 0, A9_3, 0, 0, 0, 0, 0, 0, A9_10, 0, 0, 0, 0, 0, L2],  # eq9
                [0, 0, 0, 0, 0, A10_6, 0, 0, 0, 0, A10_11, 0, L1, 0, 0, 0],  # eq10
                [0, 0, 0, 0, 0, 0, A11_7, 0, 0, 0, 0, A11_12, 0, L2, 0, 0],  # eq11
                [0, A12_2, 0, 0, 0, 0, 0, 0, A12_9, 0, 0, 0, 0, 0, L1, 0],  # eq12
                [0, 0, A13_3, 0, 0, 0, 0, 0, 0, A13_10, 0, 0, -1, 0, 0, 1],  # eq13
                [0, 0, 0, 0, 0, A14_6, 0, 0, 0, 0, A14_11, 0, 1, -1, 0, 0],  # eq14
                [0, 0, 0, 0, 0, 0, A15_7, 0, 0, 0, 0, A15_12, 0, 1, -1, 0],  # eq15
                [0, A16_2, 0, 0, 0, 0, 0, 0, A16_9, 0, 0, 0, 0, 0, 1, -1],  # eq16
            ]

            # RHS vector assembly
            b = [
                [0],
                [0],
                [0],
                [0],
                [-self.__t2],
                [-self.__t3],
                [-self.__t4],
                [-self.__t1],
                [-self.__t2 * L2 * L2 / (2 * self.__e_modulus * self.__Iyy_12)],
                [-self.__t3 * L1 * L1 / (2 * self.__e_modulus * self.__Iyy_23)],
                [-self.__t4 * L2 * L2 / (2 * self.__e_modulus * self.__Iyy_34)],
                [-self.__t1 * L1 * L1 / (2 * self.__e_modulus * self.__Iyy_41)],
                [-self.__t2 * L2 / (self.__e_modulus * self.__Iyy_12)],
                [-self.__t3 * L1 / (self.__e_modulus * self.__Iyy_23)],
                [-self.__t4 * L2 / (self.__e_modulus * self.__Iyy_34)],
                [-self.__t1 * L1 / (self.__e_modulus * self.__Iyy_41)]
            ]

            # Solve matrix problem
            # x = [Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Mz_12, Mz_23, Mz_34, Mz_41, th_ab_b, th_bc_c, th_cd_d, th_da_a]
            x = solve(A, b)

            # Transfer to class variables
            self.__F_c1_x = float(-x[1])  # F_c1_x = -Ay
            self.__F_c1_y = float(-x[0])  # F_c1_y = -Ax
            self.__F_c2_x = float(-x[3])  # F_c2_x = -By
            self.__F_c2_y = float(-x[2])  # F_c2_y = -Bx
            self.__F_c3_x = float(-x[5])  # F_c3_x = -Cy
            self.__F_c3_y = float(-x[4])  # F_c3_y = -Cx
            self.__F_c4_x = float(-x[7])  # F_c4_x = -Dy
            self.__F_c4_y = float(-x[6])  # F_c4_y = -Dx
            self.__Mz_12 = float(x[8])
            self.__Mz_23 = float(x[9])
            self.__Mz_34 = float(x[10])
            self.__Mz_41 = float(x[11])

            if ver:
                print("Ax = " + str(round(float(x[0]), 3)) + str(" N"))
                print("Ay = " + str(round(float(x[1]), 3)) + str(" N"))
                print("Bx = " + str(round(float(x[2]), 3)) + str(" N"))
                print("By = " + str(round(float(x[3]), 3)) + str(" N"))
                print("Cx = " + str(round(float(x[4]), 3)) + str(" N"))
                print("Cy = " + str(round(float(x[5]), 3)) + str(" N"))
                print("Dx = " + str(round(float(x[6]), 3)) + str(" N"))
                print("Dy = " + str(round(float(x[7]), 3)) + str(" N"))
                print("Mz_12 = " + str(round(float(x[8]), 3)) + str(" Nm"))
                print("Mz_23 = " + str(round(float(x[9]), 3)) + str(" Nm"))
                print("Mz_34 = " + str(round(float(x[10]), 3)) + str(" Nm"))
                print("Mz_41 = " + str(round(float(x[11]), 3)) + str(" Nm"))
                print("th_ab|b = " + str(round(float(x[12]) * 180 / pi, 3)) + str(" deg"))
                print("th_bc|c = " + str(round(float(x[13]) * 180 / pi, 3)) + str(" deg"))
                print("th_cd|d = " + str(round(float(x[14]) * 180 / pi, 3)) + str(" deg"))
                print("th_da|a = " + str(round(float(x[15]) * 180 / pi, 3)) + str(" deg"))

        else:
            print("Error: cptype unsupported")
            exit()

    def nvm12(self, d):
        if not 0 <= d <= self.__l2:
            print("Error: Invalid distance in nvm12()")
            exit()

        N = 0

        Mh_1 = - self.__M_c1_y
        Mh_2 = self.__M_c2_y
        Vv = (Mh_2 - Mh_1) / self.__l2

        Mv_1 = self.__Mz_12
        Mv_2 = self.__t2 + self.__Mz_23
        Vh = (Mv_2 - Mv_1) / self.__l2

        Mh = Mh_1 + Vv * d
        Mv = Mv_1 + Vh * d
        return N, Vv, Vh, Mh, Mv

    def nvm23(self, d):
        if not 0 <= d <= self.__l1:
            print("Error: Invalid distance in nvm23()")
            exit()

        N = 0

        Mh_2 = self.__M_c2_x
        Mh_3 = - self.__M_c3_x
        Vv = (Mh_3 - Mh_2) / self.__l1

        Mv_2 = self.__Mz_23
        Mv_3 = self.__t3 + self.__Mz_34
        Vh = (Mv_3 - Mv_2) / self.__l1

        Mh = Mh_2 + Vv * d
        Mv = Mv_2 + Vh * d
        return N, Vv, Vh, Mh, Mv

    def nvm34(self, d):
        if not 0 <= d <= self.__l2:
            print("Error: Invalid distance in nvm34()")
            exit()

        N = 0

        Mh_3 = self.__M_c3_y
        Mh_4 = - self.__M_c4_y
        Vv = (Mh_4 - Mh_3) / self.__l2

        Mv_3 = - self.__Mz_34
        Mv_4 = self.__t4 + self.__Mz_41
        Vh = (Mv_4 - Mv_3) / self.__l2

        Mh = Mh_3 + Vv * d
        Mv = Mv_3 + Vh * d
        return N, Vv, Vh, Mh, Mv

    def nvm41(self, d):
        if not 0 <= d <= self.__l1:
            print("Error: Invalid distance in nvm41()")
            exit()

        N = 0

        Mh_4 = - self.__M_c4_x
        Mh_1 = self.__M_c1_x
        Vv = (Mh_1 - Mh_4) / self.__l1

        Mv_4 = self.__Mz_41
        Mv_1 = self.__t1 + self.__Mz_12
        Vh = (Mv_1 - Mv_4) / self.__l1

        Mh = Mh_4 + Vv * d
        Mv = Mv_4 + Vh * d
        return N, Vv, Vh, Mh, Mv

    def tau_f(self, x, Vh, Vv):
        if not - self.__b * 0.5 <= x <= self.__b * 0.5:
            print("Error: Invalid distance x-value in tau_f()")
            exit()

        if x < 0:
            x = -x

        return Vv * (self.__b * self.__h * self.__tf / 4 - self.__h * self.__tf / 2 * x) / (self.__Ixx * self.__tf) + \
               Vh * (self.__b * self.__b * self.__tf / 8 - self.__tf / 2 * x * x) / (self.__Iyy * 2 * self.__tf)

    def tau_w(self, y, Vh, Vv):
        if not - self.__h * 0.5 <= y <= self.__h * 0.5:
            print("Error: Invalid y-value in tau_w()")
            exit()
        return Vv * (self.__h * self.__b * self.__tf / 2 + self.__h * self.__h * self.__tw
                     / 8 - self.__tw / 2 * y * y) / self.__Ixx * self.__tw

    def sigma_f(self, x, Mh, Mv):
        # Assuming constant sigma due to Mh
        # Validity: small flange thickness
        sigma = Mh * self.__h * 0.5 / self.__Ixx + Mv * x / self.__Iyy
        return sigma

    def sigma_w(self, y, Mh, Mv):
        # Assuming no sigma due to Mv
        # Validity: small distance from neutral axis
        sigma = Mh * y / self.__Ixx
        return sigma

    def localcheck(self, d, f_nvm, plot=False, print_fail=False):
        x_f = np.linspace(-self.__b / 2, self.__b / 2, 25)
        y_w = np.linspace(-self.__h / 2, self.__h / 2, 25)

        vM_f = np.zeros_like(x_f)
        vM_w = np.zeros_like(y_w)

        N, Vv, Vh, Mh, Mv = f_nvm(d)

        fail = False

        # flange
        for i in range(len(x_f)):
            # von Mises failure criterion
            vM_f[i] = sqrt(self.sigma_f(x_f[i], Mh, Mv) ** 2 + 3 * self.tau_f(x_f[i], Vh, Vv) ** 2)
            if vM_f[i] > self.__sigma_yield:
                fail = True
                if print_fail:
                    print("yield in flange")

            # compression buckling
            sigma_cr = 0.407 * pi * pi * self.__e_modulus / (12 * (1 - 0.33 * 0.33)) * (self.__tf * 2 / self.__b) ** 2
            if self.sigma_f(x_f[i], Mh, Mv) > sigma_cr:
                fail = True
                if print_fail:
                    print("flange column buckling")

        # web
        for i in range(len(y_w)):
            # von Mises failure criterion
            vM_w[i] = sqrt(self.sigma_w(y_w[i], Mh, Mv) ** 2 + 3 * self.tau_w(y_w[i], Vh, Vv) ** 2)
            if vM_f[i] > self.__sigma_yield:
                fail = True
                if print_fail:
                    print("yield in web")

            # shear buckling
            if print_fail and self.tau_w(y_w[i], Vh, Vv) > 8.2 * self.__e_modulus * (self.__tw / self.__b) ** 2:
                fail = True
                if print_fail:
                    print("web shear buckling")

        max_vM = np.max([vM_f, vM_w])

        if plot:
            plt.figure(1)
            plt.plot(x_f, vM_f / 1e6)
            plt.xlabel("x_f")
            plt.ylabel("vM_f")
            plt.figure(2)
            plt.plot(vM_w / 1e6, y_w)
            plt.xlabel("vM_w")
            plt.ylabel("y_w")
            plt.show()

        return max_vM, fail

    def globalcheck(self, prnt=False):
        l1_arr = np.linspace(0, self.__l1, 100)
        l2_arr = np.linspace(0, self.__l2, 100)

        max_vM_12 = np.zeros_like(l2_arr)
        max_vM_23 = np.zeros_like(l1_arr)
        max_vM_34 = np.zeros_like(l2_arr)
        max_vM_41 = np.zeros_like(l1_arr)

        fail = False
        for i in range(len(l2_arr)):
            max_vM_12[i], fail = self.localcheck(l2_arr[i], self.nvm12)
            if fail:
                break
            max_vM_34[i], fail = self.localcheck(l2_arr[i], self.nvm34)
            if fail:
                break

        for i in range(len(l1_arr)):
            max_vM_23[i], fail = self.localcheck(l2_arr[i], self.nvm23)
            if fail:
                break
            max_vM_41[i], fail = self.localcheck(l2_arr[i], self.nvm41)
            if fail:
                break

        globalmax = np.max([max_vM_12, max_vM_23, max_vM_34, max_vM_41])

        if prnt:
            print("von Mises - global max: " + str(round(globalmax, 2)) + " MPa")
            if globalmax > self.__sigma_yield / 1e6:
                print("Result: Yield")
            else:
                print("Result: No yield")

        return globalmax, fail

    def mass(self):
        return self.__mass


# Constants
mass_MHM = 2200  # kg
mass_motors = 110  # kg
hover_thrust_perc = 40
oei_torque = 700  # Nm

# Limit cases
safety_factor = 1.5
max_load = - mass_MHM * 9.81 * 100 / (4 * hover_thrust_perc) * safety_factor
min_load = mass_motors * 9.81 * safety_factor
max_torque = oei_torque * safety_factor
hover_thrust = mass_motors * 9.81 / 4

# Object generation
cp = Centerpiece(1.72, 1.72, 1.72, 1.72, "beams_square")

# Range of test values
b_min = 0.045
b_max = 0.055
b_arr = np.linspace(b_min, b_max, 10)

tf_min = 5e-3
tf_max = 7e-3
tf_arr = np.linspace(tf_min, tf_max, 10)

tw_min = 0.4e-3
tw_max = 1e-3
tw_arr = np.linspace(tw_min, tw_max, 1)

t_min = 2e-3
t_max = 6e-3
t_arr = np.linspace(t_min, t_max, 10)

# array to hold optimal values
# b, tf, tw, mass
opt = [0, 0, 0, 999]

# Material stats
# Al 6061-T6:
#   E_mod       = 68.9e9
#   rho         = 2700
#   sigma_yield = 260e6

# Al 2024-T4:
#   E_mod       = 73.1e9
#   rho         = 2780
#   sigma_yield = 310e6

# Al 7075-T6:
#   E_mod       = 71.7e9
#   rho         = 2810
#   sigma_yield = 448e6

# Optimal value generation
i = 0
for b in b_arr:
    print(round(b, 3))
    for t in t_arr:
        print("   " + str(round(t, 5)), end=" ")
        cp.set_beam_config(1.65, 1.65, 0.3, b, t, t, 71.7e9, 2810, 448e6)

        # OEI at full thrust
        cp.set_p(max_load, max_load, max_load, max_load)
        cp.set_t(max_torque, 0, 0, 0)

        f = False
        for theta in np.linspace(0, pi, 90):
            cp.set_arm_angle(theta, -theta, pi*0.5, pi*1.5)
            cp.solve_c()

            if cp.globalcheck()[1]:
                f = True
                print("FAIL")
                break

        if f:
            continue
        elif cp.mass() < opt[3]:
            opt = [b, t, t, cp.mass()]
            print("NEW OPT")
        else:
            print("")
        break
    i += 1
    # print(round(i / len(b_arr) * 100, 1))

print(opt)

# Results
# Al 6061-T6:
# [0.088889, 0.00544444, 0.000400000, 19.38640]
# [0.077778, 0.00444445, 0.004444445, 36.08000]

# Al 2024-T4:
# [0.083334, 0.00488889, 0.000400000, 17.15198]
# [0.077778, 0.00377778, 0.003777778, 31.57669]

# Al 7075-T6:
# [0.065556, 0.00444445, 0.000400000, 13.03258]
# [0.052222, 0.00377778, 0.003777778, 28.33646]



