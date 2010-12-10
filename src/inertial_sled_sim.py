"""
inertial_sled_sim.py

Purpose: provides a simulation of the inertial sled system with the goal of
determining the required bandwith for the short range laser distance sensor,
the the required update frequency of the realtime loop and the required motor
performance.

"""
import scipy
import scipy.integrate
import scipy.interpolate
import pylab
from force_model import ForceModel

class InertialSledSim(object):
    """
    Encapsulates a simulation of the inertial sled system. A PID controller
    is used to try and keep the sled in position with the model. 
    """

    def __init__(self,param,inertial_dt=0.1,pgain=5.0):
        # Setup force model
        self.param = param
        self.force_model = ForceModel(param)
        self.drag_func = self.force_model.get_square_model()
        # Initialize External force function
        self.force_func = zero_force_func
        # Time step of realtime system integrating the measured forces
        self.inertial_dt = inertial_dt
        self.pgain = pgain

    def normal_system_func(self,w,t):
        """
        Function governing normal system dynamics. Basic second order system with
        square law fluid dynamic drag.
        """
        v = w[0]
        x = w[1]
        m = param.sub_body_mass + param.sub_mount_mass
        dvdt = -scipy.sign(v)*self.drag_func(v)/m + self.force_func(v,t)/m
        dxdt = v
        return scipy.array([dvdt,dxdt])

    def inertial_system_func(self,x,t,vel_setpt):
        """
        Function governing the dynamics of the sled system. We are setting sled
        velocity so this is first order in position.  Currenlty doesn't
        incorporate sled/controller dynamics assumes we can control the sleds
        velocity with reasonable accuracy. Will need to add sled + controller
        dynamics later. Currently, just integration of the velocity setput.
        """
        dx = vel_setpt
        return dx

    def run_normal(self,x0,v0,T,num=100):
        """
        Runs a simulation of the system reponse to the external force function
        self.force_func. This is the normal system response - not the captive
        trajectory model.
        """
        w0 = scipy.array([v0,x0])
        t = scipy.linspace(0,T,num)
        w = scipy.integrate.odeint(self.normal_system_func,w0,t)
        x = w[:,1]
        v = w[:,0]
        return t,x,v

    def run_inertial(self,x0,v0,t_sled,x_sled,v_sled,T,num=30):
        """ 
        Run a simulation of inertial sled system. Used a proportional 
        controller with feed forward velocity term.
        """
        t_sled = t_sled - t_sled[0]
        x_sled_func = scipy.interpolate.interp1d(t_sled,x_sled)
        v_sled_func = scipy.interpolate.interp1d(t_sled,v_sled)

        t_curr = 0.0
        x_curr = x0
        v_setpt = v0

        t = []
        x = []
        v = []

        while t_curr < T:

            # Update velocity set point
            x_error = x_sled_func(t_curr) - x_curr
            v_setpt = self.pgain*x_error + v_sled_func(t_curr)

            t_next = t_curr + self.inertial_dt
            if t_next > T:
                t_next = T
            t_sim = scipy.linspace(0,t_next-t_curr,num)
            x_sim = scipy.integrate.odeint(self.inertial_system_func,x_curr,t_sim,args=(v_setpt,))

            # Add simulation segments to t and velcity lists
            N = t_sim.shape[0]
            t_sim = [val + t_curr for val in t_sim]
            t.extend(list(t_sim))
            v.extend(list(v_setpt*scipy.ones((N,))))
            x.extend(list(x_sim))

            # Update t_curr and v_curr
            t_curr = t_next
            x_curr = x_sim[-1]

       # Convert lists to arrays
        t = scipy.array(t)
        x = scipy.array(x)
        v = scipy.array(v)
        return t,x,v

    def get_perror_est(self,f0,f1,num=5):
        """
        Computes the maximum absolute positon error over a range of realtime loop update
        frequencies form f0 to f1 with the number of step given by num.
        """
        f_array = scipy.linspace(f0,f1,num)
        dt_array = 1.0/f_array
        perror_array = scipy.zeros(f_array.shape)

        for i, dt in enumerate(dt_array):

            print 'f: ', 1.0/dt

            # Run simulation for given realtime loop update interval
            self.inertial_dt = dt
            x0 = 0.0
            v0 = 0.0
            t_model, x_model, v_model= self.run_normal(x0,v0,T,num=1000)
            t_sled, x_sled, v_sled = self.run_inertial(x0,v0,t_model,x_model,v_model,T)

            # Interpolate normal trajectory onto captive trajectory time base
            x_model= scipy.reshape(x_model, t_model.shape)
            x_sled = scipy.reshape(x_sled, t_sled.shape)
            interp_func = scipy.interpolate.interp1d(t_model,x_model)
            x_model_interp = interp_func(t_sled)

            # Compute maximum absolute velocity error
            perror = scipy.absolute((x_sled - x_model_interp))
            perror_array[i] = perror.max()

        return f_array, perror_array


def zero_force_func(v,t):
    """
    Zero force function - returns zero for times and velocities.
    """
    return 0.0

def get_const_force_func(force_val):
    """
    Returns a constant force function with constant value force_val.
    """
    def const_force_func(v,t):
        return force_val
    return const_force_func

# -----------------------------------------------------------------------------
if __name__ == '__main__':

    import parameters as param

    x0 = 0.0
    v0 = 0.0
    v_step = 1.0
    T = 10.0
    num = 100

    if 1:
        # -------------------------------------------------------------------------
        # Example simulations of the inertial sled system
        # -------------------------------------------------------------------------
        inertial_f_list = [10,20,50]
        for f in inertial_f_list:

            inertial_dt = 1.0/f 
            sim = InertialSledSim(param)
            sim.inertial_dt = inertial_dt
            const_force = sim.drag_func(v_step)
            sim.force_func = get_const_force_func(const_force)
            t_model, x_model, v_model = sim.run_normal(x0,v0,T,num=num)
            t_sled, x_sled, v_sled = sim.run_inertial(x0,v0,t_model,x_model,v_model,T)

            pylab.figure()
            pylab.subplot(211)
            pylab.plot(t_model,x_model,'b')
            pylab.plot(t_sled,x_sled,'r')
            pylab.ylabel('position (m)')
            pylab.title('Inertial Sled Sim, f=%d (Hz)'%(int(f),))

            pylab.subplot(212)
            pylab.plot(t_model,v_model,'b')
            pylab.plot(t_sled,v_sled,'r')
            pylab.ylabel('velocity (m/s)')
            pylab.xlabel('time (sec)')
            pylab.savefig('iner_sled_sim_f_%d_Hz.png'%(int(f),))


    if 1:
        # -------------------------------------------------------------------------
        # Estimation of percent position error as a function of the realtime update
        # frequency.
        # -------------------------------------------------------------------------
        f0 = 10
        f1 = 500 
        num = 30 
        sim = InertialSledSim(param)
        const_force = sim.drag_func(v_step)
        sim.force_func = get_const_force_func(const_force)
        f_array, perror_array = sim.get_perror_est(f0,f1,num=num)

        pylab.figure()
        pylab.semilogy(f_array, 1000*perror_array,'o')
        pylab.grid('on')
        pylab.xlabel('realtime loop freq. (Hz)')
        pylab.ylabel('max. abs. position error (mm)') 
        pylab.title('Position Error vs Loop Freq.')
        pylab.savefig('iner_sled_sim_perror.png')
        pylab.show()


