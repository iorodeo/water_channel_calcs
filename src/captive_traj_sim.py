"""
captive_traj_sim.py

Purpose: Provides a simulation of the performance of the captive trajecory
system to help determine the bandwidth requirements for the long range laser
distance sensor and the required update frequency of system's realtime loop.
Simulations of normal (inertial) system response are also provided for
comparison. Plots example trajectories for constant force response. Plots the
percent absolute velocity error as a function of realtime loop frequency. 

"""
import scipy
import scipy.integrate
import scipy.interpolate
import pylab
from force_model import ForceModel


class CaptiveTrajSim(object):

    """
    Encapsulates a simulation of a captive trajectory system modeling a
    the first order dynamics of the sled system.

    Simulates normal and captive trajectory system for comparison. System
    update_dt can be varied to determine the effect of the realtime system
    update interval.
    """

    def __init__(self,param,captive_dt=0.5):
        # Setup force model
        self.param = param
        self.force_model = ForceModel(param)
        self.drag_func = self.force_model.get_square_model()
        # Initialize External force function
        self.force_func = zero_force_func
        # Time step of realtime system integrating the measured forces
        self.captive_dt = captive_dt

    def normal_system_func(self,v,t):
        """
        Function governing normal system dynamics, i.e., f(v,t) where

        dv/dt = f(v,t).
        """
        # Get mass and damping
        m = param.sub_body_mass + param.sub_mount_mass
        # Compute function value
        val = -scipy.sign(v)*self.drag_func(v)/m + self.force_func(v,t)/m
        return val 

    def captive_system_func(self,v,t,v_cap):
        """
        Function governing the captive trajectory system dynamics , i.e., f(v,t) where

        dv/dt = f(v,t).

        Note, In this case v is ignored and v_cap is used as the forces remain
        constant between realtime system updates. The system is simulated by
        running an odeint integration for the system evolution between every
        realtime system update.
        """
        # Get mass and damping
        m = param.sub_body_mass + param.sub_mount_mass
        # Compute function value
        val = -scipy.sign(v_cap)*self.drag_func(v_cap)/m + self.force_func(v_cap,t)/m
        return val


    def run_normal(self,v0,T,num=1000):
        """
        Runs a simulation of the system reponse to the external force function
        self.force_func. This is the normal system response - not the captive
        trajectory model.

        Inputs:

        v0  = initial velocity
        dt  = sample interval for output values
        T   = duration of simulation 
        """
        t = scipy.linspace(0,T,num)
        v = scipy.integrate.odeint(self.normal_system_func,v0,t)
        return t,v

    def run_captive(self,v0,T,num=100):
        """
        Run a simulation of the captive trajectory system response to the external 
        force function self.force_func. This is a prediction of the captive trajectory 
        system response. 
        v0 = 
        """
        t_curr = 0.0
        v_curr = v0
        t = []
        v = []

        while t_curr < T:

            # Integrate system forward to next realtime update
            t_next = t_curr + self.captive_dt
            if t_next > T:
                t_next = T
            t_sim = scipy.linspace(0,t_next-t_curr,num)
            v_sim = scipy.integrate.odeint(self.captive_system_func,v_curr,t_sim,args=(v_curr,))

            # Add simulation segments to t and velcity lists
            t_sim = [val + t_curr for val in t_sim]
            t.extend(t_sim)
            v.extend(v_sim)

            # Update t_curr and v_curr
            t_curr = t_next
            v_curr = v_sim[-1]

        # Convert lists to arrays
        t = scipy.array(t)
        v = scipy.array(v)

        return t,v


    def get_verror_est(self,f0,f1,num=5):
        """
        Get velocity error estimate as a function of the captive trajectory system
        update frequency.

        Note, the force_func must be set to something reasonable for this to work.
        A step function is a good choice.
        """
        capt_f_array = scipy.linspace(f0,f1,num)
        capt_dt_array = 1.0/capt_f_array 
        verror_array = scipy.zeros(capt_dt_array.shape)
        
        for i, capt_dt in enumerate(capt_dt_array):

            # Run simulation for given captive trajectory update interval
            self.captive_dt = capt_dt
            v0 = 0.0
            t_norm, v_norm = self.run_normal(v0,T,num=2000)
            t_capt, v_capt = self.run_captive(v0,T,num=100)

            # Interpolate normal trajectory onto captive trajectory time base
            v_norm = scipy.reshape(v_norm, t_norm.shape)
            v_capt = scipy.reshape(v_capt, t_capt.shape)
            interp_func = scipy.interpolate.interp1d(t_norm,v_norm)
            v_norm_interp = interp_func(t_capt)

            # Compute maximum absolute velocity error
            verror = scipy.absolute((v_capt - v_norm_interp))
            verror_array[i] = verror.max()

        return capt_f_array, verror_array
            

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

    v0 = 0.0
    v_step = 1.0
    T = 10.0
    num_norm = 2000
    num_capt = 200

    sim = CaptiveTrajSim(param)
    const_force = sim.drag_func(v_step)
    sim.force_func = get_const_force_func(const_force)

    # -------------------------------------------------------------------------
    # Example simulations for two different update intervals
    # -------------------------------------------------------------------------
    captive_f_list = [2.0,30.0]
    for captive_f in captive_f_list:
        sim.captive_dt = 1.0/captive_f
        t_norm, v_norm = sim.run_normal(v0,T,num=num_norm)
        t_capt, v_capt = sim.run_captive(v0,T,num=num_capt)

        pylab.figure()
        pylab.plot(t_norm,v_norm,'b')
        pylab.plot(t_capt,v_capt,'r')
        pylab.xlabel('time (sec)')
        pylab.ylabel('velocity (m/s)')
        pylab.title('Captive trajectory sim., realtime f=%d (Hz)'%(int(captive_f),))
        v_max = max([v_norm.max(), v_capt.max()])
        pylab.ylim(0,1.05*v_max)
        pylab.grid('on')
        pylab.savefig('capt_traj_sim_f_%d_Hz.png'%(int(captive_f)))

    # -------------------------------------------------------------------------
    # Estimation of percent velocity error as a function of the realtime update
    # frequency. 
    # -------------------------------------------------------------------------
    f0 = 1.5
    f1 = 100.0
    num_f = 30 
    f_vals, verror = sim.get_verror_est(f0,f1,num_f)
    verror_percent = 100.0*verror/v_step

    pylab.figure()
    pylab.semilogy(f_vals,verror_percent,'o')
    pylab.xlabel('realtime update freq. (Hz)')
    pylab.ylabel('percent velocity error')
    pylab.grid('on')
    pylab.savefig('capt_traj_sim_verror.png')
    pylab.show()

