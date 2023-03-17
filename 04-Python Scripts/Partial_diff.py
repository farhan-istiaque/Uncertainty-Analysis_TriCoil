# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 12:08:22 2022

@author: farhan
"""

"""
_____________________________________________________________
Versions:
    
python 3.9.13
numpy 1.21.5
panda 1.4.4
scipy 1.9.1
scikit-learn 1.0.2
matplotlib 3.5.2
CoolProp 6.4.1
_____________________________________________________________
"""

import math
from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI


def central_diff(var, func, thres):
    """
    This function returns the finite diff estimate of
    the partial derivative of function func at x by a threshold thres

    Inputs:
    var : float
        variables that the derivative is taken with
    func: function
        function that the derivative is taken with
    thres: float
        increment that used to evaluate the finite diff

    Outputs:
    Derivative of func at var by central diff

    """
    return (func(var+thres)-func(var-thres))/(thres*2)


def air_props_un(prop,T_a_db, un_T_db, T_a_wb, un_T_wb, P_a, un_P_a, thres=1e-5):
    '''Function to calclute the uncertainty in the air properties for given
    dry, wet bulb temperatures, and pressure

    Inputs:
    prop is the property that needs to be calculated
    T_a_db, un_T_db: air dry bul temperature and its uncertainty, K
    T_a_wb, un_T_wb: air wet bul temperature and its uncertainty, K
    P_a: Pressure (ambient), Pa
    thres: threshold for the delta h (see above docstring)

    Outputs:
    percent_xxx: percent contribution of each input to the property
    Uncer: uncertainty in the property
    '''

    # function for the wet bulb temperature
    def func_Twb(T_wb):
        # if density, need reciprocal because of calculation of specific vol.
        if prop=='Dens':
            return(1/(HAPropsSI('V','T',T_a_db,'P',P_a,'B',T_wb)))
        else:
            # for all the other properties
            return(HAPropsSI(prop,'T',T_a_db,'P',P_a,'B',T_wb))

    # function for the dry bulb temperature
    def func_Tdb(T_db):
        # if density, need reciprocal because of calculation of specific vol.
        if prop=='Dens':
            return(1/(HAPropsSI('D','T',T_db,'P',P_a,'B',T_a_wb)))
        else:
            # for all the other properties
            return(HAPropsSI(prop,'T',T_db,'P',P_a,'B',T_a_wb))

    # function for the pressure
    def func_P_a(P):
        # if density, need reciprocal because of calculation of specific vol.
        if prop=='Dens':
            return(1/(HAPropsSI('V','T',T_a_db,'P',P,'B',T_a_wb)))
        else:
            # for all the other properties
            return(HAPropsSI(prop,'T',T_a_db,'P',P,'B',T_a_wb))

    # uncertainty in property w.r.t to each input
    dp_dTb=central_diff(T_a_db, func_Tdb, thres*T_a_db)
    dp_dTw=central_diff(T_a_wb, func_Twb, thres*T_a_wb)
    dp_dP=central_diff(P_a, func_P_a, thres*P_a)

    # overall uncertainty in property
    Uncer=math.sqrt(dp_dTb**2*un_T_db**2+dp_dTw**2*un_T_wb**2+\
                    dp_dP**2*un_P_a**2)

    # percent contribution of each uncertainty
    percent_T_db=(dp_dTb**2*un_T_db**2)*100/Uncer**2
    percent_T_wb=(dp_dTw**2*un_T_wb**2)*100/Uncer**2
    percent_P_a=(dp_dP**2*un_P_a**2)*100/Uncer**2

    return(Uncer, percent_T_db, percent_T_wb, percent_P_a)

def air_props_un_humid(prop,T_a_db, un_T_db, T_a_wb, un_T_wb, P_a, un_P_a, thres=1e-5):
    '''Function to calclute the uncertainty in the air properties for given
    dry, wet bulb temperatures, and pressure

    Inputs:
    prop is the property that needs to be calculated
    T_a_db, un_T_db: air dry bul temperature and its uncertainty, K
    T_a_wb, un_T_wb: air wet bul temperature and its uncertainty, K
    P_a: Pressure (ambient), Pa
    thres: threshold for the delta h (see above docstring)

    Outputs:
    percent_xxx: percent contribution of each input to the property
    Uncer: uncertainty in the property
    '''

    # function for the wet bulb temperature
    def func_Twb(T_wb):
        # if density, need reciprocal because of calculation of specific vol.
        if prop=='Dens':
            return(1/(HAPropsSI('V','T',T_a_db,'P',P_a,'B',T_wb)))
        else:
            # for all the other properties
            return(HAPropsSI(prop,'T',T_a_db,'P',P_a,'B',T_wb))

    # function for the dry bulb temperature
    def func_Tdb(T_db):
        # if density, need reciprocal because of calculation of specific vol.
        if prop=='Dens':
            return(1/(HAPropsSI('D','T',T_db,'P',P_a,'B',T_a_wb)))
        else:
            # for all the other properties
            return(HAPropsSI(prop,'T',T_db,'P',P_a,'B',T_a_wb))

    # function for the pressure
    def func_P_a(P):
        # if density, need reciprocal because of calculation of specific vol.
        if prop=='Dens':
            return(1/(HAPropsSI('V','T',T_a_db,'P',P,'B',T_a_wb)))
        else:
            # for all the other properties
            return(HAPropsSI(prop,'T',T_a_db,'P',P,'B',T_a_wb))

    # uncertainty in property w.r.t to each input
    dp_dTb=central_diff(T_a_db, func_Tdb, thres*T_a_db)
    dp_dTw=central_diff(T_a_wb, func_Twb, thres*T_a_wb)
    dp_dP=central_diff(P_a, func_P_a, thres*P_a)

    # overall uncertainty in property
    Uncer=math.sqrt(dp_dTb**2*un_T_db**2+dp_dTw**2*un_T_wb**2+\
                    dp_dP**2*un_P_a**2)

    # percent contribution of each uncertainty
    percent_T_db=(dp_dTb**2*un_T_db**2)*100/Uncer**2
    percent_T_wb=(dp_dTw**2*un_T_wb**2)*100/Uncer**2
    percent_P_a=(dp_dP**2*un_P_a**2)*100/Uncer**2

    return(Uncer, percent_T_db, percent_T_wb, percent_P_a)


def ref_props_un(prop, T_r, un_T_r, P_r, un_P_r, fluid, thres=1e-5):
    '''Function to calclute the uncertainty in the refrigerant properties for 
    given temperature and pressure

    Inputs:
    prop is the property that needs to be calculated
    T_r, un_T_r: refrigerant inlet temperature and its uncertainty, K
    P_r, un_P_r: refrigerant inlet pressure and its uncertainty, Pa
    fluid: refrigerant
    thres: threshold for the delta h (see above docstring)

    Outputs:
    Uncer: uncertainty in the property
    '''
    # function for the temperature
    def func_T_r(T):
        return(PropsSI(prop, 'T', T, 'P', P_r, fluid))

    # function for the pressure
    def func_P_r(P):
        return(PropsSI(prop, 'T', T_r, 'P', P, fluid))
    
    # uncertainty in property w.r.t to each input
    dh_dT=central_diff(T_r, func_T_r, thres*T_r)
    dh_dP=central_diff(P_r, func_P_r, thres*P_r)
    
    # overall uncertainty in property
    Uncer=math.sqrt(dh_dT**2*un_T_r**2+dh_dP**2*un_P_r**2)
    percent_T_r=(dh_dT**2*un_T_r**2)*100/Uncer**2
    percent_P_r=(dh_dP**2*un_P_r**2)*100/Uncer**2
    
    
    return(Uncer, percent_T_r, percent_P_r)

"""
TESTING THE FUNCITONS
"""
if __name__=='__main__':
    T_a_db=20+273.15
    T_a_wb=14+273.15
    prop='H'
    T_db_un=0.7
    T_wb_un=0.5
    P_a=101325
    P_a_un=30.94
    
    A=air_props_un(prop,T_a_db,T_db_un,T_a_wb,T_wb_un,P_a,P_a_un)

