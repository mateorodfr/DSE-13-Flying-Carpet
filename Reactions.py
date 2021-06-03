from math import sin, cos


class Reactions:
    def __init__(self, l_a1, l_a2, l_a3, l_a4):
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
        """

        # ----------------------
        # Initializing variables
        # ----------------------

        self.l_a1 = l_a1
        self.l_a2 = l_a2
        self.l_a3 = l_a3
        self.l_a4 = l_a4

        # They are empty but may be used later
        # If they are used in calculations before they are set, an error is returned
        # This is nice since it prevents the thrusts from being used before they are specified explicitly

        # Initializing engine thrust variables
        self.p1 = None  # p1 is (from entrance perspective) the front left engine
        self.p2 = None  # p2 is (from entrance perspective) the back left engine
        self.p3 = None  # p3 is (from entrance perspective) the back right engine
        self.p4 = None  # p4 is (from entrance perspective) the front right engine

        # Initializing engine torque variables
        self.t1 = None  # t1 is (from entrance perspective) the front left engine
        self.t2 = None  # t2 is (from entrance perspective) the front left engine
        self.t3 = None  # t3 is (from entrance perspective) the front left engine
        self.t4 = None  # t4 is (from entrance perspective) the front left engine

        # Initializing arm angles
        self.th1 = None  # th1 is (from entrance perspective) the front left arm
        self.th2 = None  # th2 is (from entrance perspective) the front left arm
        self.th3 = None  # th3 is (from entrance perspective) the front left arm
        self.th4 = None  # th4 is (from entrance perspective) the front left arm

    def set_p(self, p1, p2, p3, p4):
        """
        Sets the force produced at each arm tip
        This includes the weight of the engine
        Downward positive (positive z-direction)
        """
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

    def set_t(self, t1, t2, t3, t4):
        """
        Sets the torque of each engine
        Positive according to right hand rule in z-direction
        """
        self.t1 = t1
        self.t2 = t2
        self.t3 = t3
        self.t4 = t4

    def set_arm_angles(self, th1, th2, th3, th4):
        """
        Sets the angles of the engine arms
        Angles are relative to x-axis
        Positive according to right hand rule along z-axis
        """
        self.th1 = th1
        self.th2 = th2
        self.th3 = th3
        self.th4 = th4

    def get_h1(self):
        """
        Returns the forces and moments acting on hinge h1
        """
        # Assume thrust always acts along z-axis : rotation of motor assumed negligible
        # F_h1_x and F_h1_y are therefore neglected
        # Contribution of thrust force to M_h1_z is therefore also neglected
        F_h1_z = self.p1
        M_h1_x = self.p1 * self.l_a1 * sin(self.th1)
        M_h1_y = - self.p1 * self.l_a1 * cos(self.th1)
        M_h1_z = self.t1
        return F_h1_z, M_h1_x, M_h1_y, M_h1_z

    def get_h2(self):
        """
        Returns the forces and moments acting on hinge h2
        """
        # Assume thrust always acts along z-axis : rotation of motor assumed negligible
        # F_h2_x and F_h2_y are therefore also zero
        # Contribution of thrust force to M_h2_z is therefore also neglected
        F_h2_z = self.p2
        M_h2_x = self.p2 * self.l_a2 * sin(self.th2)
        M_h2_y = - self.p2 * self.l_a2 * cos(self.th2)
        M_h2_z = self.t2
        return F_h2_z, M_h2_x, M_h2_y, M_h2_z

    def get_h3(self):
        """
        Returns the forces and moments acting on hinge h3
        """
        # Assume thrust always acts along z-axis : rotation of motor assumed negligible
        # F_h3_x and F_h3_y are therefore also zero
        # Contribution of thrust force to M_h3_z is therefore also neglected
        F_h3_z = self.p3
        M_h3_x = self.p3 * self.l_a3 * sin(self.th3)
        M_h3_y = - self.p3 * self.l_a3 * cos(self.th3)
        M_h3_z = self.t3
        return F_h3_z, M_h3_x, M_h3_y, M_h3_z

    def get_h4(self):
        """
        Returns the forces and moments acting on hinge h4
        """
        # Assume thrust always acts along z-axis : rotation of motor assumed negligible
        # F_h4_x and F_h4_y are therefore also zero
        # Contribution of thrust force to M_h4_z is therefore also neglected
        F_h4_z = self.p4
        M_h4_x = self.p4 * self.l_a4 * sin(self.th4)
        M_h4_y = - self.p4 * self.l_a4 * cos(self.th4)
        M_h4_z = self.t4
        return F_h4_z, M_h4_x, M_h4_y, M_h4_z
