#!/usr/bin/env python3
# ============================================================================
# File:     cable_2P1.py
# Author:   Erik Johannes Husom
# Created:  2019-10-10
# ----------------------------------------------------------------------------
# Description:
# 
# Third mandatory exercise in IN5270, University of Oslo, fall 2019.
#
# ----------------------------------------------------------------------------
# DESCRIPTION OF PROBLEM
#
# Compute deflection of a cable with sine functions. We have a hanging cable
# with tension, and the cable has a deflection w(x) which is governed by:
#
# (1) Tw''(x) = l(x),
#
# - L: Length of cable
# - T: Tension on cable
# - w(x): Deflection of cable
# - l(x): Vertical load per unit length
# 
# Cable fixed at x = 0 and x = L. Boundary conditions w(0) = w(L) = 0. Deflection 
# is positive upwards and l is positive when it acts downwards.
#
# Assuming l(x) = const, the solution is symmetric around x = L/2. For a
# function w(x) that is symmetric around a point x_0, we have that
#
# (2) w(x_0 - h) = w(x_0 + h),
#
# which means that
#
# (3) lim_{h->0} (w(x_0+h) - w(x_0 - h))/(2h) = 0.
#
# We can therefore halve the domain, since it is symmetric. That limits the
# problem to find w(x) in [0, L/2], with boundary conditions w(0) = 0 and
# w'(L/2) = 0.
#
# Scaling of variables:
#
# x_ = x/(L/2)      (setting x = x_ in code for easier notation)
# u = w/w_c         (where w_c is a characteristic size of w)
#
# By putting this into equation (1) we get
#
# (4) (4Tw_c)/L^2 * u''(x_) = l = const. 
#
# We set |u''(x_)|= 1, and we get w_c = 0.25lL^2/T, and the scaled problem is
#
# (5) u'' = 1, x_ \in (0,1), u(0) = 0, u'(1) = 0.
#
# ----------------------------------------------------------------------------
# EXACT SOLUTION FOR DEFLECTION u
#
#
# ============================================================================



