"""
force_calcs.py

Purpose: Estimate of fluid drag on the submarine model and mount -  will be
used to help select the force sensor required or if more than one sensor is
needed to span the full 0-3m/s range. Note, this is only the passive forces not
that due to the actuator. However, if the sub can really go 3m/s per second
than the actuator must be capable of producing at least this much force.

Note, I estimate the forces acting on the body and mount separately. This
should be a bit of an over estimate. Forces in N are converted to kg for  force
sensor sizing.

"""
import scipy

def get_body_force(density,drag_coef,diam,vel):
    """
    Computes the drag force acting an object with a circular
    cross section. 
    """
    return 0.5*density*drag_coef*scipy.pi*((0.5*diam)**2)*vel**2


def get_mount_force(density,drag_coef,chord,length,vel):
    """
    Computes the drag forces acting on an air foil section with
    given drag coeff, chord and length. 
    """

    return 0.5*density*drag_coef*chord*length*vel**2

# -----------------------------------------------------------------------------
if __name__ == "__main__":

    import pylab
    import parameters as param

    # Number of velocity points
    num_vel_pts = 100

    # Get drag forces acting on model body and mount
    min_vel, max_vel = param.sub_velocity_range
    vel_array = scipy.linspace(min_vel,max_vel,num_vel_pts)
    body_force_array = get_body_force(
            param.water_density,
            param.sub_body_drag_coef, 
            param.sub_body_diameter,
            vel_array
            )
    mount_force_array = get_mount_force(
            param.water_density,
            param.sub_mount_drag_coef,
            param.sub_mount_chord,
            param.sub_mount_length,
            vel_array
            )
    total_force_array = body_force_array + mount_force_array

    # Convert force arrays to kg for sensor sizing
    body_force_array_kg = body_force_array/param.grav_accel
    mount_force_array_kg = mount_force_array/param.grav_accel
    total_force_array_kg = total_force_array/param.grav_accel
    
    # Plot forces
    pylab.figure(1)
    h0 = pylab.plot(vel_array, body_force_array, 'b', label= 'body force')
    h1 = pylab.plot(vel_array, mount_force_array, 'g', label = 'mount force')
    h2 = pylab.plot(vel_array, total_force_array, 'r', label = 'total force')
    pylab.title('Sub model forces (N) vs velocity (m/s)')
    pylab.xlabel('velocity (m/s)')
    pylab.ylabel('force (N)')
    pylab.grid('on')
    pylab.legend()
    pylab.savefig('forces_N.png')

    pylab.figure(2)
    h0 = pylab.plot(vel_array, body_force_array_kg, 'b', label = 'body force')
    h1 = pylab.plot(vel_array, mount_force_array_kg, 'g', label = 'mount force')
    h2 = pylab.plot(vel_array, total_force_array_kg, 'r', label = 'total force')
    pylab.title('Sub model body forces (kg) vs velocity (m/s)')
    pylab.xlabel('velocity (m/s)')
    pylab.ylabel('force (kg)')
    pylab.grid('on')
    pylab.legend()
    pylab.savefig('force_kg.png')

    pylab.show()


