"""
parameters.py

Purpose: provides a single location for storing parameters relating to
calculations for the Caltech 40m water channel project.

Note, some numbers are necessarily estimates and should be clearly marked as
such. In general I have tried to select all estimated parameters such that they
over estimate the forces acting on the model. 

Author: Will Dickson
"""

# Constants
INCH2M = 0.0254
LB2KG = 0.45359237 
water_density = 1000.0
grav_accel = 9.81

# Submarine model body parameters
sub_body_diameter = 6.0*INCH2M
sub_body_mass = 30.0*LB2KG  # Estimate
sub_body_drag_coef= 1.0     # Estimate

# Submarine model mount (hydrofoil) parameters
sub_mount_chord = 5.0*INCH2M
sub_mount_length = 14.5*INCH2M
sub_mount_max_width = 1.0*INCH2M  # Estimate
sub_mount_mass = 10*LB2KG         # Estimate
sub_mount_drag_coef = 0.2         # Estimate

# Submaring motion parameters
sub_velocity_range  = (0.0, 3.0)


