
"""
moment_calcs.py

Purpose: estimate of fluid dynamic moments on the submarine model and mount -
will be used to help size the air bearing and sensor moment and side load
requirements for this application.  

Note, only calculates the passive moments not those due to the actuator We will
need discuss the actuator thrust with John and Robert. I assume the sub itseld
is symetrical and doesn't generate any moment in normal operation.  Thus moment
should be due to the hydrofoil connecting the sub to the sled. I assume the CM of
the sub is located along the subs center line. 

"""
import scipy


def get_mount_moment(density, drag_coef, chord, r0, r1, vel):
    """
    Computes moment due to the hydrofoil conecting the sub model to the sled.
    The estimate comes from integrating the moments rxf along the mount
    hydrofoil. 

    r0 = start of mount from model CM
    r1 = end of mount form model CM
    """
    return 0.25*density*drag_coef*chord*(vel**2)*(r1**2 - r0**2)


if __name__ == '__main__':

    import pylab
    import parameters as param

    # Number of velocity points
    num_vel_pts = 100

    # Get drag forces acting on model body and mount
    min_vel, max_vel = param.sub_velocity_range
    vel_array = scipy.linspace(min_vel,max_vel,num_vel_pts)
    r0 = 0.5*param.sub_body_diameter
    r1 = param.sub_mount_length + r0
    moment_array = get_mount_moment(
            param.water_density, 
            param.sub_mount_drag_coef, 
            param.sub_mount_chord,
            r0, 
            r1,
            vel_array
            )

    # Convert to kg m to help with sizing of force sensor 
    moment_array_kg = moment_array/param.grav_accel

    # Plot results
    pylab.figure(1)
    pylab.plot(vel_array, moment_array)
    pylab.xlabel('velocity (m/s)')
    pylab.ylabel('moment (N m)')
    pylab.grid('on')
    pylab.savefig('moment_N_m.png')

    pylab.figure(2)
    pylab.plot(vel_array, moment_array_kg)
    pylab.xlabel('velocity (m/s)')
    pylab.ylabel('moment (kg m)')
    pylab.grid('on')
    pylab.savefig('moment_kg_m.png')
    pylab.show()


