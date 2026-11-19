"""Integrity visualizer for schema-3 trajectories: a GIF film and a static PNG panel. Not a paper figure.

The film shows saved configurations only (every k-th saved frame, no interpolation) in the lab frame, styled after
Fig. 2: body gray, leading tip black, centre-of-mass trail red. The panel shows s(t), the COM path, the |s| histogram,
the railway MSD of paper bead j=25 against a passive reference, and the mean-square rotation. All numbers in the panel
come from spinchain.analysis, which reads only the HDF5 file.
"""
from __future__ import annotations
import json,math,os,h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from .analysis import analyze
from .io import strict_json

TITLE='integrity visualization, not a paper figure'
NO_SPIRAL='no spiral in this window'  # written when max |s| < 1 over the saved frames; a reporting note, not a regime label

def select_frames(n,max_frames=200):
    """Indices of saved frames used in the film: every k-th frame with k=ceil(n/max_frames). No interpolation."""
    return np.arange(0,n,max(1,math.ceil(n/max_frames)))

def load(path):
    with h5py.File(path,'r') as f:
        a=f.attrs
        if 'origin_convention' not in a:raise ValueError(f'{path}: origin_convention attribute missing; refusing to draw an unlabelled frame')
        if 'positions' not in f:raise ValueError(f'{path}: no positions dataset; rerun headless without --no-positions')
        d={k:np.asarray(f[k]) for k in ('step','time','time_tau','spiral','com','positions')}
        d['pars']=json.loads(a['parameters_json']);d['attrs']={k:(v.item() if isinstance(v,np.generic) else v) for k,v in a.items() if k!='parameters_json'}
    return d

def com_reset_check(d):
    """Evidence that the stored COM was never reset: it equals the mean of the stored lab-frame positions, and no
    frame after the first has |COM| ~ 0 (the old recenter wrote ~1e-16 there)."""
    com=d['com'];later=d['step']>0
    eq=bool(np.array_equal(com,d['positions'].mean(axis=1)));zero=int(np.sum(np.linalg.norm(com[later],axis=1)<1e-9))
    return {"com_equals_mean_of_positions":eq,"frames_after_start_with_abs_com_below_1e-9":zero,"com_reset_detected":(not eq) or zero>0,
        "com_first":com[0].tolist(),"com_last":com[-1].tolist(),"com_displacement":float(np.linalg.norm(com[-1]-com[0]))}

def _label(p):return rf"Pe={p['pe']:g}   $\xi_P/L$={p['xi_over_l']:g}   N={p['n_bonds']}   seed={p['seed']}"

def film(d,out,fps=12,max_frames=200,title=TITLE):
    p=d['pars'];pos=d['positions'];com=d['com'];idx=select_frames(len(pos),max_frames);lead=pos.shape[1]-1  # code bead n_beads-1 leads
    pts=pos[idx].reshape(-1,2);lo=pts.min(0);hi=pts.max(0);c=(lo+hi)/2;half=0.55*float(np.max(hi-lo))+2
    fig,ax=plt.subplots(figsize=(6,6.6),dpi=90);fig.suptitle(title,fontsize=11)
    ax.set_xlim(c[0]-half,c[0]+half);ax.set_ylim(c[1]-half,c[1]+half);ax.set_aspect('equal')
    ax.set_xlabel(r'x / $r_0$ (lab frame)');ax.set_ylabel(r'y / $r_0$')
    trail,=ax.plot([],[],color='red',lw=1.2,label='centre of mass')
    body,=ax.plot([],[],color='0.6',lw=2.5,marker='o',ms=2.5,mfc='0.6',mec='none',label='filament')
    tip,=ax.plot([],[],'o',color='black',ms=6,label='leading tip')
    ax.legend(loc='lower right',fontsize=8,frameon=False)
    txt=ax.text(0.02,0.98,'',transform=ax.transAxes,va='top',fontsize=9)
    def draw(k):
        i=idx[k];x=pos[i];body.set_data(x[:,0],x[:,1]);tip.set_data([x[lead,0]],[x[lead,1]]);trail.set_data(com[:i+1,0],com[:i+1,1])
        txt.set_text(_label(p)+'\n'+rf"t = {d['time_tau'][i]:.3e} $\tau$   s = {d['spiral'][i]:+.2f}"+f"\nsaved frame {i} (step {int(d['step'][i]):,}), film frame {k+1}/{len(idx)}")
        return body,tip,trail,txt
    if np.max(np.abs(d['spiral']))<1:ax.text(0.02,0.02,NO_SPIRAL+' (max |s| < 1)',transform=ax.transAxes,fontsize=9)
    FuncAnimation(fig,draw,frames=len(idx)).save(out,writer=PillowWriter(fps=fps));plt.close(fig)
    return {"film_frames":int(len(idx)),"film_stride_saved_frames":int(idx[1]-idx[0]) if len(idx)>1 else 1,"fps":fps}

def panel(d,z,out,check,title=TITLE):
    p=d['pars'];pos=d['positions'];com=d['com'];lead=pos.shape[1]-1
    fig,ax=plt.subplots(2,3,figsize=(16,9.5),dpi=100);fig.suptitle(f"{title}  |  "+_label(p),fontsize=12)
    a=ax[0,0];a.plot(d['time_tau'],d['spiral'],color='k',lw=1);a.axhline(0,color='0.8',lw=0.8)
    a.ticklabel_format(axis='x',style='sci',scilimits=(0,0));a.set_xlabel(r'$t/\tau$');a.set_ylabel('spiral number s');a.set_title('s(t), eq. (12)')
    a=ax[0,1];idx=select_frames(len(pos),4)
    for k,i in enumerate(idx):
        g=str(0.75-0.45*k/max(1,len(idx)-1));a.plot(pos[i,:,0],pos[i,:,1],color=g,lw=1.5);a.plot(pos[i,lead,0],pos[i,lead,1],'o',color='k',ms=4)
    a.plot(com[:,0],com[:,1],color='red',lw=1.2);a.set_aspect('equal',adjustable='datalim')
    a.set_xlabel(r'x / $r_0$ (lab frame)');a.set_ylabel(r'y / $r_0$');a.set_title('COM path (red), saved snapshots (gray, tip black)')
    a=ax[0,2];h=z['spiral']['abs_s_histogram'];e=np.array(h['edges']);a.bar(e[:-1],h['density'],width=h['bin_width'],align='edge',color='0.5',edgecolor='k',lw=0.5)
    a.set_xlabel('|s|');a.set_ylabel('p(|s|)');a.set_title('|s| histogram (bin 0.25)')
    b2=z['spiral']['kurtosis_beta2_symmetrized'];a.text(0.97,0.95,(r'$\beta_2$ (symmetrized) = %.3g'%b2) if b2 is not None else r'$\beta_2$ undefined',transform=a.transAxes,ha='right',va='top',fontsize=9)
    a=ax[1,0];r=z.get('railway');v=p['f_p']*p['r0']/p['gamma']
    if r:
        lt=np.array(r['lag_time']);cols=plt.cm.viridis(np.linspace(0,0.85,len(r['msd_j_plus_i'])))
        for (i,m),cl in zip(r['msd_j_plus_i'].items(),cols):
            a.loglog(lt[1:],np.array(m)[1:],color=cl,label=f'bead j+{i}')
            tr=int(i)*p['r0']/v if v>0 else np.inf
            if tr<=lt[-1]:a.axvline(tr,color=cl,ls=':',lw=0.8)
        pr=r.get('passive_reference')
        if pr:a.loglog(np.array(pr['lag_time'])[1:],np.array(pr['self_msd_bead_j'])[1:],color='k',lw=2,label='passive bead j (Pe=0)')
        else:a.text(0.03,0.05,'no passive reference supplied',transform=a.transAxes,fontsize=9)
        miss=[i for i in r['msd_j_plus_i'] if v<=0 or int(i)*p['r0']/v>lt[-1]]
        if miss:a.text(0.03,0.95,'transit i r0/v beyond max lag for i='+','.join(miss)+':\nrailway incomplete in this window',transform=a.transAxes,va='top',fontsize=8)
        a.legend(fontsize=8,loc='lower right',frameon=False)
    a.set_xlabel('lag t (simulation time units)');a.set_ylabel(r'MSD / $r_0^2$');a.set_title(f"railway MSD, paper bead j={r['paper_bead_j'] if r else 25} (Fig. 6 style); dotted: i r0/v")
    a=ax[1,1];rt=z['rotation'];lt=np.array(rt['lag_time']);msr=np.array(rt['msr']);ok=(lt>0)&(msr>0);a.loglog(lt[ok],msr[ok],color='b',label='MSR')
    fit=rt['dr_fit_unit_slope']
    if fit:
        tt=np.geomspace(fit['t_min'],fit['t_max'],20);a.loglog(tt,2*fit['D']*tt,'r--',label='unit-slope fit (inspect)')
        a.text(0.03,0.95,f"mechanical unit-slope fit on t={fit['t_min']:.3g}..{fit['t_max']:.3g}\nD_r={fit['D']:.3g} (inspect before use)",transform=a.transAxes,va='top',fontsize=8)
    else:a.text(0.03,0.95,'no unit-slope region: D_r not reported',transform=a.transAxes,va='top',fontsize=9)
    a.legend(fontsize=8,loc='lower right',frameon=False);a.set_xlabel('lag t (simulation time units)');a.set_ylabel(r'MSR / rad$^2$');a.set_title('mean-square rotation of r_e (Fig. 8a style)')
    a=ax[1,2];a.axis('off');w=z['window']
    lines=[z['replication_status'],'',f"steps {int(d['step'][-1]):,}, saved frames {len(d['step'])}, sample every {d['attrs'].get('sample_every')} steps",
        f"window {w['start_tau']:.3g}..{w['end_tau']:.3g} tau ({w['reason']})",f"max |s| = {float(np.max(np.abs(d['spiral']))):.3g}"+(f"  -> {NO_SPIRAL}" if np.max(np.abs(d['spiral']))<1 else ''),
        f"COM reset detected: {check['com_reset_detected']}",f"  COM == mean(positions): {check['com_equals_mean_of_positions']}",
        f"  frames with |COM|<1e-9 after start: {check['frames_after_start_with_abs_com_below_1e-9']}",f"  COM displacement {check['com_displacement']:.3g} r0",
        '',f"policy {p['nonbonded_policy']}: hypothesis, unconfirmed",f"dt {p['dt']:g}: unpublished (code guard 0.2 gamma/k_S)",
        "origin: lab frame, never recentered","leading tip: code bead N = paper bead 1",f"tau = {p['tau']:g} time units; lags in time units"]
    a.text(0,1,'\n'.join(lines),va='top',fontsize=9.5,family='monospace',transform=a.transAxes)
    fig.tight_layout(rect=(0,0,1,0.96));fig.savefig(out);plt.close(fig)

def visualize(path,out_dir='output/viz',passive_path=None,fps=12,max_frames=200,window_tau=(0.0,None),window_reason='visualization: whole run, transient included',title=TITLE):
    d=load(path);z=analyze(path,include_structure=True,window_tau=window_tau,window_reason=window_reason,passive_path=passive_path)
    os.makedirs(out_dir,exist_ok=True);base=os.path.join(out_dir,os.path.splitext(os.path.basename(path))[0]);check=com_reset_check(d)
    info=film(d,base+'.gif',fps,max_frames,title);panel(d,z,base+'.png',check,title)
    s=np.abs(d['spiral']);out={"trajectory":str(path),"passive_reference":None if passive_path is None else str(passive_path),"gif":base+'.gif',"png":base+'.png',
        "steps":int(d['step'][-1]),"saved_frames":int(len(d['step'])),"sample_every":d['attrs'].get('sample_every'),**info,
        "max_abs_s":float(s.max()),"duration_tau":float(d['time_tau'][-1]),"com_check":check,"d_r_fit":z['rotation']['dr_fit_unit_slope'],
        "replication_status":z['replication_status']}
    strict_json(base+'.json',out);return out

def plateau_film(path,out,first,last,fps=12,max_frames=200,title=TITLE,margin=3.0,follow=False,trail_from_first=False):
    """Film of saved frames first..last only (inclusive), boxed on the chain over those frames, not on the whole COM
    trail. Lab-frame coordinates; every frame is a saved configuration (stride k=ceil(n/max_frames), no interpolation).
    follow=True: a fixed-size box re-centred on each frame's COM, for plateaus over which the COM wanders far.
    trail_from_first=True: the red COM trail starts at frame `first` instead of the start of the run."""
    d=load(path);p=d['pars'];pos=d['positions'];com=d['com'];lead=pos.shape[1]-1
    ct=p['length']*p['gamma']*(p['n_bonds']+1)/p['length']/p['f_p'];tc=d['time']/ct
    idx=np.arange(first,last+1)[select_frames(last-first+1,max_frames)]
    pts=pos[idx].reshape(-1,2);lo=pts.min(0);hi=pts.max(0);c=(lo+hi)/2;half=0.5*float(np.max(hi-lo))+margin
    if follow:half=float(max(np.max(np.abs(pos[i]-com[i])) for i in idx))+margin
    fig,ax=plt.subplots(figsize=(6,6.6),dpi=90);fig.suptitle(title,fontsize=11)
    ax.set_xlim(c[0]-half,c[0]+half);ax.set_ylim(c[1]-half,c[1]+half);ax.set_aspect('equal')
    ax.set_xlabel(r'x / $r_0$ (lab frame, zoomed on the chain)');ax.set_ylabel(r'y / $r_0$')
    trail,=ax.plot([],[],color='red',lw=1.2,label='centre of mass')
    body,=ax.plot([],[],color='0.6',lw=2.5,marker='o',ms=3,mfc='0.6',mec='none',label='filament')
    tip,=ax.plot([],[],'o',color='black',ms=7,label='leading tip')
    ax.legend(loc='lower right',fontsize=8,frameon=False);txt=ax.text(0.02,0.98,'',transform=ax.transAxes,va='top',fontsize=9)
    def draw(k):
        i=idx[k];x=pos[i];body.set_data(x[:,0],x[:,1]);tip.set_data([x[lead,0]],[x[lead,1]]);t0=first if trail_from_first else 0
        trail.set_data(com[t0:i+1,0],com[t0:i+1,1])
        if follow:ax.set_xlim(com[i,0]-half,com[i,0]+half);ax.set_ylim(com[i,1]-half,com[i,1]+half)
        txt.set_text(_label(p)+'\n'+f"t = {tc[i]:.3f} contours   s = {d['spiral'][i]:+.2f}"+f"\nsaved frame {i} (step {int(d['step'][i]):,}), film frame {k+1}/{len(idx)}")
        return body,tip,trail,txt
    FuncAnimation(fig,draw,frames=len(idx)).save(out,writer=PillowWriter(fps=fps));plt.close(fig)
    return {"gif":str(out),"saved_frames":[int(first),int(last)],"film_frames":int(len(idx)),"box_half_width_r0":half,"follow":bool(follow),
        "t_contours":[float(tc[first]),float(tc[last])]}
