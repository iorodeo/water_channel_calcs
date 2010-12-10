"""
force_model.py

Purpose: provides a models of the forces and moments acting on the submarine
model and mount hydrofoil. Both the standard square law model and linearized
models force models are given.  The linearized force models are specified by
selecting a desired operating point. Note, this module uses functions in the
force_calcs module.

"""
import scipy
import force_calcs

class ForceModel(object):

    def __init__(self, param):
        self.param = param

    def get_linear_coef(self,v0):
        """
        Returns linearized model coefficients: F(v0) and F'(v0).
        """
        # Get reference areas for sub body and mount
        sub_body_area = scipy.pi*(0.5*self.param.sub_body_diameter)**2
        sub_mount_area = self.param.sub_mount_length*self.param.sub_mount_chord

        # Get sum of product of force coefficients and reference areas
        coef_x_area = sub_body_area*self.param.sub_body_drag_coef
        coef_x_area += sub_mount_area*self.param.sub_mount_drag_coef
        F_v0 = 0.5*self.param.water_density*coef_x_area*v0**2
        dF_v0 = self.param.water_density*coef_x_area*v0
        return F_v0, dF_v0

    def get_square_model(self):
        """
        Returns a force function for the sub body + mount based on a square law
        force model. 
        """
        def square_model(v):
            f_body = force_calcs.get_body_force(
                    self.param.water_density, 
                    self.param.sub_body_drag_coef, 
                    self.param.sub_body_diameter,
                    v
                    )
            f_mount = force_calcs.get_mount_force( 
                    self.param.water_density, 
                    self.param.sub_mount_drag_coef, 
                    self.param.sub_mount_chord, 
                    self.param.sub_mount_length, 
                    v
                    )

            return f_body+f_mount

        return square_model


    def get_linear_model(self,v0):
        """
        Returns a force function for the sub body + mount based on the
        linearized force model about the operating point v0.

        Model of form F(V) =  F(v0) + F'(V0)*v0
        """

        F_v0, dF_v0 = self.get_linear_coeff(v0)

        def linear_model(v):
            return F_v0 + dF_v0*(v-v0)

        return linear_model


# ------------------------------------------------------------------------------
if __name__ == '__main__':

    """
    Test linearized model at several test points to ensure that is implemented
    correctly. 
    """
    import pylab
    import parameters as param

    num_test_pts = 5
    num_vel_pts = 100
    dv0 = 0.5
    min_vel, max_vel = param.sub_velocity_range
    v0_test_pts = scipy.linspace(0.1*max_vel, max_vel,5)
    force_model  =  ForceModel(param)

    # Plot results for square law model
    square_model = force_model.get_square_model()
    v = scipy.linspace(min_vel,max_vel, num_vel_pts)
    f_square = square_model(v)
    pylab.plot(v,f_square)

    for v0 in v0_test_pts:
        linear_model = force_model.get_linear_model(v0)

        v = scipy.linspace(v0-dv0,v0+dv0,num_vel_pts)
        f_linear = linear_model(v)
        pylab.plot(v,f_linear)

    pylab.xlabel('velocty (m/s)')
    pylab.ylabel('force (N)')
    pylab.grid('on')
    pylab.savefig('linear_model_test.png')
    pylab.show()

