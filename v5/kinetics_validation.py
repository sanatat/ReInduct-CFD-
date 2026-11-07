"""Measured NiO product-layer growth-rate validation with composition holdout.

Hidayat et al. (2009), Metallurgical and Materials Transactions B 40,
Table III: initial Ni interface layer growth in micrometres per second.
This surface-growth observable is NOT equal to the bulk volumetric rate
in the reactor code; no unsupported unit conversion is performed.
https://doi.org/10.1007/s11663-008-9212-0
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Exact original Table III, measured reaction rates and reported 1-sigma SD.
# Each row (temperature_C, H2 fraction of H2-N2 gas, rate_um_s, std_um_s)
ROWS=[
 (500.,1.,.61,.09), (600.,1.,2.08,.05), (675.,1.,4.04,.12),
 (800.,1.,2.80,.10), (900.,1.,3.82,.10), (1000.,1.,5.25,.40),
 (675.,.75,2.66,.06), (800.,.75,1.94,.13),
 (900.,.75,2.35,.07), (1000.,.75,3.54,.10),
 (600.,.53,.72,.04), (675.,.53,1.97,.21), (800.,.53,1.46,.12),
 (900.,.53,1.77,.12), (1000.,.53,2.45,.15),
 (675.,.26,1.15,.03), (800.,.26,.70,.09),
 (900.,.26,.89,.13), (1000.,.26,1.34,.05)]

def design(T_C, x_h2, mode='segmented'):
    T_C=np.asarray(T_C,float);x_h2=np.asarray(x_h2,float)
    z=1000*(1./(T_C+273.15)-1./1000.)
    if mode=='arrhenius':return np.column_stack([np.ones_like(z),z,np.log(x_h2)])
    if mode=='segmented':
        low=(T_C<=725).astype(float);high=1.-low
        return np.column_stack([low,low*z,high,high*z,np.log(x_h2)])
    raise ValueError(mode)

def fit(T,y,sd,p,mode,order_fixed=None):
    X=design(T,p,mode)
    logy=np.log(y); sigma=np.maximum(sd/y,.025)
    if order_fixed is not None:
        X=X[:,:-1]
        logy=logy-order_fixed*np.log(p)
    a=X/sigma[:,None];b=logy/sigma
    coef=np.linalg.lstsq(a,b,rcond=None)[0]
    if np.linalg.matrix_rank(a)<X.shape[1]:raise ValueError('underdetermined fit')
    return coef

def predict(coef,T,p,mode,order_fixed=None):
    X=design(T,p,mode)
    if order_fixed is not None:
        return np.exp(np.clip(X[:,:-1]@coef+order_fixed*np.log(p),-100,100))
    return np.exp(np.clip(X@coef,-100,100))

def assess(out='v5/results/kinetics'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    a=np.asarray(ROWS,float);T,p,rate,sd=a.T
    train=p>=.53; test=p==.26
    metrics={}
    predictions={}
    for label,mode,order in [('arrhenius','arrhenius',None),('segmented_free_order','segmented',None),('segmented_first_order','segmented',1.)]:
        coef=fit(T[train],rate[train],sd[train],p[train],mode,order)
        pred=predict(coef,T,p,mode,order)
        testerrors=pred[test]-rate[test]
        metrics[label]={'fit_coefficients':coef.tolist(), 'hydrogen_reaction_order':('fitted' if order is None else 1.),
          'heldout_H2_fraction':.26,'heldout_points':int(test.sum()),
          'heldout_MAE_um_s':float(np.mean(np.abs(testerrors))),
          'heldout_RMSE_um_s':float(np.sqrt(np.mean(testerrors**2))),
          'heldout_normalized_RMSE_vs_reported_sd':float(np.sqrt(np.mean((testerrors/sd[test])**2))),
          'train_weighted_RMSE':float(np.sqrt(np.mean(((pred[train]-rate[train])/sd[train])**2)))}
        predictions[label]=pred
    # Measurement-uncertainty-only bootstrap; tiny dataset, no systematic covariance accounted for.
    rng=np.random.default_rng(2606); vals=[]
    for k in range(1200):
        sim=np.maximum(.01,rate[train]+rng.normal(0.,sd[train]))
        cc=fit(T[train],sim,sd[train],p[train],'segmented',1.)
        vals.append(predict(cc,T[test],p[test],'segmented',1.))
    vals=np.asarray(vals)
    lo,hi=np.quantile(vals,[.025,.975],axis=0)
    metrics['segmented_first_order']['predictive_parametric_MC95_band_um_s']={'temperatures_C':T[test].tolist(),
                  'lo':lo.tolist(),'hi':hi.tolist(), 'samples':1200}
    metrics['source']='Hidayat et al 2009 Table III DOI 10.1007/s11663-008-9212-0'
    metrics['experiment']='Dense synthetic NiO reduced in H2-N2 at 1 atm; initial 0-5 s layer-growth rate'
    metrics['training']='Concentrations 53%, 75%, 100% H2, all tabulated temperatures'
    metrics['heldout']='26% H2 condition, four temperatures unseen during fitting'
    metrics['caveats']=['Empirical segmented fit jumps across unobserved 725C regime boundary; first-order model was evaluated after identifying pressure holdout error (exploratory comparison)',
       'Prediction uncertainty is only from tabulated SD; not systematic/model structural uncertainty',
       'Do not translate um/s of metal product-layer growth into NiO bed mol/m3/s without geometry',
       'Independent composition holdout inside one experimental paper; NOT validated reactor conversion']
    (out/'validation.json').write_text(json.dumps(metrics,indent=2)+'\n')
    with (out/'comparison.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['temperature_C','H2_fraction','measured_growth_um_s','sd_um_s','split',
                                    'single_arrhenius_pred_um_s','two_regime_first_order_pred_um_s'])
        for i in range(len(T)):
            w.writerow([T[i],p[i],rate[i],sd[i],'holdout' if test[i] else 'fit',
                        predictions['arrhenius'][i],predictions['segmented_first_order'][i]])
    fig,axs=plt.subplots(1,2,figsize=(11,4.2))
    for ph in (1.,.75,.53):
        mask=p==ph
        axs[0].errorbar(T[mask],rate[mask],yerr=sd[mask],fmt='o',label=f'{int(ph*100)}% H2')
    axs[0].set_title('Training: literature initial Ni-layer growth');axs[0].set_xlabel('Temperature [°C]')
    axs[0].set_ylabel('Measured growth [µm/s]');axs[0].legend(fontsize=8)
    axs[1].errorbar(T[test],rate[test],yerr=sd[test],fmt='ko',label='Unseen 26% H2 measured')
    for label,style in [('arrhenius','--'),('segmented_first_order','-'),('segmented_free_order',':')]:
        axs[1].plot(T[test],predictions[label][test],style,marker='s',label=label)
    axs[1].fill_between(T[test],lo,hi,alpha=.16,label='segmented MC 95% (measurement SD only)')
    axs[1].set_title('Holdout: 26% H2 at 675–1000 °C')
    axs[1].set_xlabel('Temperature [°C]');axs[1].set_ylabel('Layer-growth rate [µm/s]')
    axs[1].legend(fontsize=7)
    for ax in axs:ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'kinetics_holdout.png',dpi=165);plt.close(fig)
    print(json.dumps(metrics,indent=2))
    return metrics

if __name__=='__main__':assess()
