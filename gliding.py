###
# gliding.py
# functions and classes to simulate 2d gliding mechanics
# by: sam hocking
# updated: 5/9/2022
###

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
from scipy.integrate import odeint
from scipy.optimize import brentq

g = 9.81
m_per_ft = 0.3048
deg_per_radian = 180/np.pi

def lin_interp(df, field, xlookup, x):
    '''
    linear interpolation function
    df: reference dataframe
    field: string formatted dataframe field name (dependent variable)
    lookup: string formatted dataframe field name (independent variable)
    x: value to interoplate on
    '''
    filtered = copy.deepcopy(df.iloc[(df[xlookup]-x).abs().argsort()[:2]])
    filtered.sort_values(xlookup, inplace=True)
    alpha_arr = np.array(filtered[xlookup])
    field_arr = np.array(filtered[field])
    return np.interp(x, alpha_arr, field_arr)

def Re_interp(Re_arr, df_arr, field, alpha, Re):
    '''
    airfoil coefficient interpolation across array of polar dataframes given an input Reynolds number
    Re_arr: array of Reynolds number values
    df_arr: array of airfoil polar dataframes corresponding to the Re_arr values
    field: string formatted dataframe field name (dependent variable)
    lookup: string formatted dataframe field name (independent variable)
    Re: Reynolds number to interoplate on
    '''
    sort_idx = np.abs((Re_arr-Re)).argsort()[:2]
    filtered_Re_arr = Re_arr[sort_idx]
    df1 = df_arr[sort_idx[0]]
    df2 = df_arr[sort_idx[1]]
    val_arr = np.array([lin_interp(df1, field, 'alpha', alpha), lin_interp(df2, field, 'alpha', alpha)])
    return np.interp(Re, filtered_Re_arr, val_arr)

class airfoil():
    '''
    airfoil property class
    imports list of polar files (cleaned exports from XFOIL) for varying Reynolds number
    polar_files: list of file names
    Re_list: list of corresponding Reynolds number values
    name: airfoil name
    '''
    def __init__(self, polar_files, Re_list, name):
        self.polar_files = polar_files
        self.Re_list = Re_list
        self.name = name
        polar_df_arr = []
        for i in range(len(self.Re_list)):
            temp_df = pd.read_csv(self.polar_files[i], header=0)
            temp_df.sort_values('alpha', inplace=True)
            temp_df['CL/CD'] = temp_df['CL'] / temp_df['CD']
            polar_df_arr.append(temp_df)
        self.polar_df_arr = polar_df_arr
        self.cm_roots = self.cm_root_find_loop(self.polar_df_arr)
    
    def cm_root_bracket(self, df):
        '''
        scan through polar dataframe and find angle of attack (alpha) values that bracket a c_m root
        '''
        pairs_arr = []
        alpha_arr = np.array(df['alpha'])
        cm_arr = np.array(df['CM'])
        for i in range(len(alpha_arr)-1):
            if cm_arr[i]*cm_arr[i+1] < 0:
                pairs_arr.append([alpha_arr[i],alpha_arr[i+1]])
        return pairs_arr

    def cm_root_find(self, df):
        '''
        use the root brackets and use linear interpolation to approximate c_m roots
        '''
        brackets = self.cm_root_bracket(df)
        roots = []
        def root_find_f(x):
            return lin_interp(df, 'CM', 'alpha', x)
        for pair in brackets:
            root = brentq(root_find_f, pair[0], pair[1])
            roots.append(root)
        return np.array(roots)

    def cm_root_find_loop(self, df_arr):
        '''
        iterate through each polar dataframe to find the c_m roots
        '''
        roots_arr = []
        for df in df_arr:
            roots_arr.append(self.cm_root_find(df))
        return roots_arr
    
    def plot(self, x_field, y_field, size=(6,6)):
        '''
        flexible polar plotting function
        '''
        fig, ax = plt.subplots(figsize=size)
        for i in range(len(self.Re_list)):
            df = self.polar_df_arr[i]
            plt.plot(df[x_field], df[y_field], label=f'Re={self.Re_list[i]}')
        ax.legend()
        ax.set_xlabel(x_field)
        ax.set_ylabel(y_field)
        plt.title(f'{y_field} vs. {x_field} [{self.name}]')
        plt.show()

    def plot_cm_fixed_pts(self, size=(6,6)):
        '''
        more specific plotting function to plot c_m vs. alpha with the c_m roots (fixed points)
        '''
        fig, ax = plt.subplots(figsize=(6,6))
        for i in range(len(self.Re_list)):
            df = self.polar_df_arr[i]
            plt.plot(df['alpha'], df['CM'], label=f'Re={self.Re_list[i]}')
        for i in range(len(self.Re_list)):
            for j in range(len(self.cm_roots[i])):
                plt.axvline(x=self.cm_roots[i][j], linestyle='--')
        ax.legend()
        ax.set_xlabel('alpha')
        ax.set_ylabel('CM')
        plt.title(f'CM vs. alpha [{self.name}]')
        plt.show()

class atmosphere():
    '''
    import the standard SI engineering atmosphere
    alt: geometric altitude (above sea level) in meters
    temp: temp in K
    press: pressure in N/m^2
    dens: density in kg/m^2
    speed_of_sound: speed of sound in m/s
    visc: dynamic viscosity of air in kg/(m*s)
    alt_ft: geometric altitude (above sea level) in feet
    '''
    def __init__(self, atmos_file):
        atmos_df = pd.read_csv(atmos_file, header=0, usecols = ['alt', 'temp', 'press', 'dens', 'a', 'visc'])
        atmos_df.rename(columns={'a': 'speed_of_sound'}, inplace=True)
        atmos_df['alt'] = atmos_df['alt'] * 1000
        atmos_df['visc'] = atmos_df['visc'] * 1e-6
        atmos_df['alt_ft'] = atmos_df['alt'] / m_per_ft
        self.df = atmos_df

class airframe():
    '''
    airframe class
    takes dimension arguments to define basic aircraft properties
    c: chord length
    cg_offest: offsetting vector from center of gravity located at [0,0] to the aerodynamic center (approximated by quarter chord)
    wingspan in m
    wing_area in m^2
    m: mass in kg 
    '''
    def __init__(self, c, cg_offset, wingspan, wing_area, m):
        self.c = c
        self.c_qtr = c/4
        self.cg_offset = cg_offset
        self.wingspan = wingspan
        self.wing_area = wing_area
        self.m = m
        self.cg = np.array([0,0])
        self.ac = self.cg + self.cg_offset
        self.cg_ac_len = np.linalg.norm(self.ac - self.cg)
        self.aspect_ratio = self.wingspan**2 / self.wing_area

class simulate():
    '''
    parent simulation class
    atmos_class: initialized atmosphere object
    airfoil_class: initialized airfoil object
    airframe_class: initialized airframe class
    vect0: initial conditions vector [v1, v2, v3, v4, v5, v6]
      v1: initial x position of center of gravity
      v2: initial y position of center of gravity
      v3: initial angle between center of gravity and the aerodynamic center and horizontal (inclination of the ac)
      v4: initial x velocity
      v5: initial y velocity
      v6: initial angular velocity of the aerodynamic center
    '''
    def __init__(self, atmos_class, airfoil_class, airframe_class, vect0):
        self.atmos = atmos_class
        self.airfoil = airfoil_class
        self.airframe = airframe_class
        self.vect0 = vect0

    def integrate(self, t0, tf, dt):
        '''
        basic integration function using odeint
        t0: initial time
        tf: final time
        dt: sampling time step
        '''
        tarr = np.arange(t0,tf,dt)
        solution = odeint(self.dvect, self.vect0, tarr)
        varr = solution[:,3]**2 + solution[:,4]**2
        self.results = np.column_stack((tarr, solution, varr))

class freefall(simulate):
    '''
    simulation object for freefall behavior
    '''
    def __init__(self, atmos_class, airfoil_class, airframe_class, vect0):
        super().__init__(atmos_class, airfoil_class, airframe_class, vect0)
        self.type='freefall'

    def dvect(self, y, t):
        '''
        given y vector, returns y' given gravity
            v1: x position of center of gravity
            v2: y position of center of gravity
            v3: angle between center of gravity and the aerodynamic center and horizontal (inclination of the ac)
            v4: x velocity
            v5: y velocity
            v6: angular velocity of the aerodynamic center        
        '''
        v1, v2, v3, v4, v5, v6 = y
        dvect = [v4, v5, v6, 0, -g, 0]
        return dvect

class stable_glide(simulate):
    '''
    simulation object for stable gliding behavior
    '''
    def __init__(self, atmos_class, airfoil_class, airframe_class, vect0):
        super().__init__(atmos_class, airfoil_class, airframe_class, vect0)
        self.type='stable glide'

    def dvect(self, y, t):
        '''
        given y vector, returns y' given gravity, lift, and drag
        assumes alpha is constant at initial alpha (perfect, constant trim)
            v1: x position of center of gravity
            v2: y position of center of gravity
            v3: angle between center of gravity and the aerodynamic center and horizontal (inclination of the ac)
            v4: x velocity
            v5: y velocity
            v6: angular velocity of the aerodynamic center        
        '''
        v1, v2, v3, v4, v5, v6 = y
        theta = np.arctan(v5 / v4)
        ac_x = v1 + self.airframe.cg_ac_len*np.cos(v3)
        ac_y = v2 + self.airframe.cg_ac_len*np.sin(v3)
        alpha = v3
        vel_sq = v4**2 + v5**2
        rho = lin_interp(self.atmos.df, 'dens', 'alt', ac_y)
        mu = lin_interp(self.atmos.df, 'visc', 'alt', ac_y)
        Re = rho * np.sqrt(v4**2 + v5**2) * self.airframe.c_qtr / mu
        cl = Re_interp(self.airfoil.Re_list, self.airfoil.polar_df_arr, 'CL', alpha*180/(np.pi), Re)
        cd = Re_interp(self.airfoil.Re_list, self.airfoil.polar_df_arr, 'CD', alpha*180/(np.pi), Re)
        L = 1/2 * rho * vel_sq * self.airframe.c * self.airframe.wingspan * cl
        D = 1/2 * rho * vel_sq * self.airframe.c * self.airframe.wingspan * cd
        dvect = [v4, v5, v6, L/self.airframe.m * np.cos(np.pi/2 + theta) + D/self.airframe.m * np.cos(np.pi + theta), -g + L/self.airframe.m * np.sin(np.pi/2 + theta) + D/self.airframe.m * np.sin(np.pi + theta), 0]
        return dvect

class dynamic_glide(simulate):
    '''
    simulation object for dynamic gliding behavior
    '''
    def __init__(self, atmos_class, airfoil_class, airframe_class, vect0):
        super().__init__(atmos_class, airfoil_class, airframe_class, vect0)
        self.type='dynamic glide'

    def dvect(self, y, t):
        '''
        given y vector, returns y' given gravity, lift, drag, and a moment about the center of gravity
        assumes the aerodynamic center is the quarter chord
        Mcg (moment about the center of gravity) = Mac + l*L*cos(alpha) - l*D*sin(alpha)
            v1: x position of center of gravity
            v2: y position of center of gravity
            v3: angle between center of gravity and the aerodynamic center and horizontal (inclination of the ac)
            v4: x velocity
            v5: y velocity
            v6: angular velocity of the aerodynamic center        
        '''
        v1, v2, v3, v4, v5, v6 = y
        theta = np.arctan(v5 / v4)
        ac_x = v1 + self.airframe.cg_ac_len*np.cos(v3)
        ac_y = v2 + self.airframe.cg_ac_len*np.sin(v3)
        alpha = theta - v3
        vel_sq = v4**2 + v5**2
        rho = lin_interp(self.atmos.df, 'dens', 'alt', ac_y)
        mu = lin_interp(self.atmos.df, 'visc', 'alt', ac_y)
        Re = rho * np.sqrt(v4**2 + v5**2) * self.airframe.c_qtr / mu
        cl = Re_interp(self.airfoil.Re_list, self.airfoil.polar_df_arr, 'CL', alpha*180/(np.pi), Re)
        cd = Re_interp(self.airfoil.Re_list, self.airfoil.polar_df_arr, 'CD', alpha*180/(np.pi), Re)
        cm = Re_interp(self.airfoil.Re_list, self.airfoil.polar_df_arr, 'CM', alpha*180/(np.pi), Re)
        L = 1/2 * rho * vel_sq * self.airframe.c * self.airframe.wingspan * cl
        D = 1/2 * rho * vel_sq * self.airframe.c * self.airframe.wingspan * cd
        Mac = 1/2 * rho * vel_sq * self.airframe.c * self.airframe.wingspan * self.airframe.c * cm
        Mcg = (Mac + self.airframe.cg_ac_len*L*np.cos(alpha) - self.airframe.cg_ac_len*D*np.sin(alpha))        
        dvect = [v4, v5, v6, L/self.airframe.m * np.cos(np.pi/2 + theta) + D/self.airframe.m * np.cos(np.pi + theta), -g + L/self.airframe.m * np.sin(np.pi/2 + theta) + D/self.airframe.m * np.sin(np.pi + theta), Mcg/(self.airframe.m*self.airframe.cg_ac_len**2)]
        return dvect