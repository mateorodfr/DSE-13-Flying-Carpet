from math import sin, cos


class Reactions:
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
            TODO: "beams cross" : 4 beams in a cross pattern
            TODO: "truss_coarse" : truss with low number of members
            TODO: "truss_fine" : truss with high number of members
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
            self.__moi12 = None
            self.__moi23 = None
            self.__moi34 = None
            self.__moi41 = None
            self.__e_modulus = None

    def set_p(self, p1, p2, p3, p4):
        """
        Sets the force produced at each arm tip
        This includes the weight of the engine
        Downward positive (positive z-direction)
        """
        self.__p1 = p1
        self.__p2 = p2
        self.__p3 = p3
        self.__p4 = p4

    def set_t(self, t1, t2, t3, t4):
        """
        Sets the torque of each engine
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

    def get_h1(self):
        """
        Returns the forces and moments acting on hinge h1
        """
        # Assume thrust always acts along z-axis : rotation of motor assumed negligible
        # F_h1_x and F_h1_y are therefore neglected
        # Contribution of thrust force to M_h1_z is therefore also neglected
        F_h1_z = self.__p1
        M_h1_x = self.__p1 * self.__l_a1 * sin(self.__th1)
        M_h1_y = - self.__p1 * self.__l_a1 * cos(self.__th1)
        M_h1_z = self.__t1
        return F_h1_z, M_h1_x, M_h1_y, M_h1_z

    def get_h2(self):
        """
        Returns the forces and moments acting on hinge h2
        """
        # Assume thrust always acts along z-axis : rotation of motor assumed negligible
        # F_h2_x and F_h2_y are therefore also zero
        # Contribution of thrust force to M_h2_z is therefore also neglected
        F_h2_z = self.__p2
        M_h2_x = self.__p2 * self.__l_a2 * sin(self.__th2)
        M_h2_y = - self.__p2 * self.__l_a2 * cos(self.__th2)
        M_h2_z = self.__t2
        return F_h2_z, M_h2_x, M_h2_y, M_h2_z

    def get_h3(self):
        """
        Returns the forces and moments acting on hinge h3
        """
        # Assume thrust always acts along z-axis : rotation of motor assumed negligible
        # F_h3_x and F_h3_y are therefore also zero
        # Contribution of thrust force to M_h3_z is therefore also neglected
        F_h3_z = self.__p3
        M_h3_x = self.__p3 * self.__l_a3 * sin(self.__th3)
        M_h3_y = - self.__p3 * self.__l_a3 * cos(self.__th3)
        M_h3_z = self.__t3
        return F_h3_z, M_h3_x, M_h3_y, M_h3_z

    def get_h4(self):
        """
        Returns the forces and moments acting on hinge h4
        """
        # Assume thrust always acts along z-axis : rotation of motor assumed negligible
        # F_h4_x and F_h4_y are therefore also zero
        # Contribution of thrust force to M_h4_z is therefore also neglected
        F_h4_z = self.__p4
        M_h4_x = self.__p4 * self.__l_a4 * sin(self.__th4)
        M_h4_y = - self.__p4 * self.__l_a4 * cos(self.__th4)
        M_h4_z = self.__t4
        return F_h4_z, M_h4_x, M_h4_y, M_h4_z

    def set_beam_config(self, l1, l2, moi12, moi23, moi34, moi41, e_modulus):
        """
        float l1 : y-distance between beams 12 and 34
        float l2 : x-distance between beams 23 and 41
        float moi12 : Moment of inertia of beam 12
        float moi23 : Moment of inertia of beam 23
        float moi34 : Moment of inertia of beam 34
        float moi41 : Moment of inertia of beam 41
        float e_modulus : E-modulus of material (same material for all beams)
        """
        if self.__cptype == "beams_square":
            self.__moi12 = moi12
            self.__moi23 = moi23
            self.__moi34 = moi34
            self.__moi41 = moi41
            self.__e_modulus = e_modulus
        else:
            print("Error: set_beam_config called while cptype not a beam type")
            exit()

    def get_c1(self):
        """
        Returns forces and moments acting on connector c1
        """

        if self.__cptype == "beams_square":
            # TODO: Assemble matrix containing equations of statically indeterminate problem
            # TODO: Solve matrix using sp.solve
            # TODO: Return F_c1_x and F_c1_y reactions at connectors


