"""
freq_response.py

Computes the expected frequency response of the system as a function of
operating point velocity. Plots the amplitude and phase of the system as a
function of the frequency - i.e., produces a bode plot of the system response.
Plots the time constant and 3dB frequencies as a function of operating point
velocity.

System:

    m*dv/dt = -b*v + f(t)

    where m is the system mass and b the linearized damping coefficient. 
    the linearized damping coefficient is given by:

    density*(cd_body*area_body + cd_mount*area_mount)*v

Note, as the drag forces acting on the system are a square of the velocity
the frequency response of the system will depend upon the velocity. 
"""
import scipy
import math


def get_freq_response(f0,f1,m,b,n=5000):
    """
    Get frequency response for system. Returns gain, phase and frequency.   

    Inputs:
    f0 = starting frequency
    f1 = stopping frequency
    m  = mass of system
    b  = damping coefficient

    Outputs:
    mag_db = gain (output/input)
    phase  = phase shift in degrees
    f      = array of frequencies
    """

    def transfer_func(s,m,b):
        return 1.0/(m*s + b)

    f = scipy.linspace(f0,f1,n)
    x = 2.0*scipy.pi*f*1j
    y = transfer_func(x,m,b)
    mag = scipy.sqrt(y.real**2 + y.imag**2)
    phase = scipy.arctan2(y.imag, y.real)
    phase = scipy.rad2deg(phase)
    mag_db = 20.0*scipy.log10(mag)
    return mag_db, phase, f


def get_time_const(v0,param):
    """
    Calculates the systems time constant as a function of the operating
    velocity. 
    """
    # Setup force model and predict gain and phase
    force_model  =  ForceModel(param)
    m = param.sub_body_mass + param.sub_mount_mass
    dummy, b = force_model.get_linear_coef(v0)
    return m/b

def get_f_3dB(v0,param):
    """
    Calculates the systems cut off frequency as a function of the operating 
    velocity.
    """
    time_const = get_time_const(v0,param)
    f_3dB = 1.0/(2.0*scipy.pi*time_const)
    return f_3dB

if __name__ == '__main__':

    """
    Get plots of the systems frequency response. Note, this is a nonlinear
    system. However, we can linearize the system about an operating point
    (velocity) and get the frequency response for the linearized system. This
    is done for a range of velocities from 5% of the max velocity to max
    velocity. 
    
    In addition the time constant and cut off frequency of the system is
    calculated as a function of the operating velocity.
    """
    import pylab
    import parameters as param
    from force_model import ForceModel 

    # --------------------------------------------------------------------------
    # Part 1 Compute example bode plots for system frequency response as a 
    # function of the operatiing point. Based on linearized force model.
    # --------------------------------------------------------------------------

    # Test frequency range - f0 must be low enough that there is no phase shift
    f0 = 0.001
    f1 = 10.0

    # Test velocities
    num_test_v0 = 2  
    min_vel, max_vel = param.sub_velocity_range
    test_v0 = scipy.linspace(0.05*max_vel, max_vel,num_test_v0)

    # Setup force model and predict gain and phase
    force_model  =  ForceModel(param)
    m = param.sub_body_mass + param.sub_mount_mass

    for i, v0 in enumerate(test_v0):

        # Compute the frequency response of the system
        dummy, b = force_model.get_linear_coef(v0)
        gain, phase, f = get_freq_response(f0,f1,m,b,n=10000)

        # Normalize gains - note, input and output don't have same units
        gain = gain - gain[0]

        # Plot results
        pylab.figure(1) 
        pylab.subplot(211)
        pylab.semilogx(f,gain,label='v0 = %1.2f m/s'%(v0,))
        pylab.ylim(gain.min(),1)

        pylab.subplot(212)
        pylab.semilogx(f,phase,label='v0 = %1.2f m/s'%(v0,))
        pylab.ylim(phase.min()-5,5)

    pylab.subplot(211)
    pylab.grid('on')
    pylab.ylabel('gain (dB)')
    pylab.title('Frequency Response')
    pylab.subplot(212)
    pylab.xlabel('frequency (Hz)')
    pylab.ylabel('phase (deg)')
    pylab.grid('on')
    pylab.legend()
    pylab.savefig('freq_response.png')

    # -------------------------------------------------------------------------
    # Part 2: compute the system time constant and cutoff frequency as a 
    # function of the operating point velocity
    # -------------------------------------------------------------------------
    v = scipy.linspace(0.05*max_vel, max_vel)
    tc = get_time_const(v,param)
    f_3dB = get_f_3dB(v,param)

    pylab.figure(2)
    pylab.plot(v,tc)
    pylab.title('Time Constant vs operating point')
    pylab.xlabel('operating point (m/s)')
    pylab.ylabel('time constant (s)')
    pylab.grid('on')
    pylab.savefig('time_const.png')
    pylab.figure(3)
    pylab.plot(v,f_3dB)
    pylab.title('Cutoff frequency vs operating point')
    pylab.xlabel('operating point (m/s)')
    pylab.ylabel('f cutoff (Hz)')
    pylab.grid('on')
    pylab.savefig('f_cutoff.png')
    pylab.show()


