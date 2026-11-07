"""Within-refresh nonlinear feedback between Joule heat and solid temperature.

The EM solver can be axisymmetric while the thermal reactor remains 1D.
Only thermal-conductivity feedback is closed at an EM refresh: gas transfer,
chemistry, geometry and radiation are lagged over the time step.  This is NOT
fully coupled MHD/plasma CFD.
"""
import numpy as np


def solve_em_thermal_feedback(solid_energy, solid_rhs_without_em, dt, heat_capacity,
                              temperature_old, em_solver, *, tolerance=1.e-7,
                              max_iter=24, relaxation=0.85):
    if not 0<relaxation<=1: raise ValueError('relaxation must be in (0,1]')
    q,report=em_solver(temperature_old)
    for count in range(1,max_iter+1):
        t_trial=(solid_energy+dt*(solid_rhs_without_em+q))/heat_capacity
        if not np.all(np.isfinite(t_trial)) or np.any(t_trial<=0):
            raise FloatingPointError('bad solid temperature in EM feedback')
        q_updated,rep_updated=em_solver(t_trial)
        error=float(np.max(np.abs(q_updated-q))/max(1.,np.max(np.abs(q_updated))))
        if error<tolerance:
            return q_updated,{**rep_updated,'feedback_iterations':count,
                              'feedback_relative_power_density_error':error,
                              'feedback_temperature_change_K':float(np.max(np.abs(t_trial-temperature_old)))}
        q=relaxation*q_updated+(1.-relaxation)*q
        report=rep_updated
    raise RuntimeError(f'EM/thermal feedback did not converge in {max_iter} iterations')
