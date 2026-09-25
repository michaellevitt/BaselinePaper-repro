"""
Figure S1, v2 (2026-09-19) — level-space companion to Figure 1, with the PRIMARY age-specific
baseline overlaid.

Per-population ln of the ESP-2013 age-standardised death rate, indexed to the fitted 2019 level,
on a common vertical axis.  The vertical gap between observed (black) and a baseline over
2020-2025 is that model's excess; 2020-2023 is shaded.

Age-specific trends became primary on 19 Sep 2026, so the three country-level curves are
relabelled as all-ages trends, and a fourth curve is added in orange: the age-standardised
aggregate of the twenty per-band primary STTa baselines,

    ln( sum_a ESP_a * exp(alpha_a) * clip(exp(b_a t + g_a t^2)) / sum_a ESP_a )

indexed to the same 2019 level.  This is not an approximation of the primary model, it is the
primary model's own baseline expressed in standardised space, so the gap between the black and
the orange curve is the excess the paper now reports.

Reads output/option1_country_params.csv and output/agespecific_band_params_v1.csv.
Writes output/fig_option1_logquad_montage_v2.png .
"""
import os, sys, csv, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from excess_anchor_window import load, cell_for
OUTD=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"output")
FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
ESP={"0":1000,"1-4":4000,"5-9":5500,"10-14":5500,"15-19":5500,"20-24":6000,"25-29":6000,"30-34":6500,"35-39":7000,"40-44":7000,
     "45-49":7000,"50-54":7000,"55-59":6500,"60-64":6000,"65-69":5500,"70-74":5000,"75-79":4000,"80-84":2500,"85-89":1500,"90+":1000}
ESPsum=sum(ESP.values()); W=5
RELIABLE={"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,"POL":2008,"SVK":2006}
COBS="#111"; CTT="#c0392b"; CST="#159090"; CSP="#6a3d9a"
data,name=load()
def DP(l,y,a): return cell_for(data,l,y,a,"T")
def has_base(l):
    se=max(2003+(W-1),RELIABLE.get(l,2003)+(W-1)); return all(DP(l,y,a) is not None for y in range(se-(W-1),2020) for a in FINE)
PC=[r[0] for r in csv.reader(open(os.path.join(OUTD,"slope_of_slopes_CI.csv"))) if r and r[0]!="country" and not r[0].startswith(("WEIGHTED","UNWEIGHTED","SD"))]
locs=sorted([l for l in sorted(PC) if has_base(l)], key=lambda l: name.get(l,l))
def esp(l,y):
    if any(DP(l,y,a) is None for a in FINE): return None
    return sum(ESP[a]*DP(l,y,a)[0]/DP(l,y,a)[1] for a in FINE)/ESPsum
def _f(x):
    try: return float(x)
    except (TypeError,ValueError): return None
par={r["country"]:r for r in csv.DictReader(open(os.path.join(OUTD,"option1_country_params.csv")))}
CBW="#e67e22"
BP={}
for r in csv.DictReader(open(os.path.join(OUTD,"agespecific_band_params_v1.csv"))):
    BP.setdefault(r["country"],{})[r["band"]]=(_f(r["islope2019_shrunk"]), _f(r["sos_shrunk"]))
CLIP=(0.4,2.5)
def agg_primary(l,start,a19):
    """Standardised aggregate of the 20 per-band primary STTa baselines, indexed to the 2019 fit."""
    if l not in BP: return None
    yrs=np.arange(start,2026); out=[]
    alph={}
    fy=np.array([y for y in range(start,2020)],float)
    for a in FINE:
        r=np.array([DP(l,int(y),a)[0]/DP(l,int(y),a)[1] for y in fy])
        alph[a]=float(np.polyfit(fy-2019.0,np.log(np.maximum(r,1e-9)),2)[2])
    for y in yrs:
        t_=y-2019.0; s=0.0
        for a in FINE:
            g19,sos=BP[l].get(a,(None,None))
            if g19 is None or sos is None: return None
            b,g=g19/100.0, sos/200.0
            m=min(max(np.exp(b*t_+g*t_*t_),CLIP[0]),CLIP[1])
            s+=ESP[a]*np.exp(alph[a])*m
        out.append(np.log(s/ESPsum)-a19)
    return yrs, np.array(out)
NR,NC=8,5
D={}; allv=[]
for l in locs:
    start=max(2003,RELIABLE.get(l,2003))
    obs_y=[y for y in range(start,2026) if esp(l,y) is not None]; obsln=[np.log(esp(l,y)) for y in obs_y]
    fit_y=[y for y in range(start,2020)]
    c=np.polyfit(np.array(fit_y,float)-2019.0,[np.log(esp(l,y)) for y in fit_y],2)   # uniform log-quad -> fitted 2019 level
    a19=c[2]
    oi={y:np.log(esp(l,y))-a19 for y in obs_y}                                        # observed, indexed
    p=par.get(l,{}); bo=_f(p.get("beta_i")); go=_f(p.get("gamma_i")); bs=_f(p.get("beta_shrunk")); gs=_f(p.get("gamma_shrunk"))
    xf=np.arange(start,2026); tf=xf-2019.0
    tta=bo*tf+go*tf*tf if bo is not None else None                                    # TTa own (indexed)
    stta=bs*tf+gs*tf*tf if bs is not None else None                                   # STTa shrunk (indexed)
    bas=_f(p.get("beta_anch_shrunk")); gas=_f(p.get("gamma_anch_shrunk")); sttap=None
    if bas is not None and gas is not None:                                            # STTa+ = slope-anchored shrunk trend (parabola)
        xp=np.arange(2019,2026); tp=xp-2019.0; sttap=(xp, bas*tp+gas*tp*tp)
    prim=agg_primary(l,start,a19)
    D[l]=(start,obs_y,oi,xf,tta,stta,sttap,prim)
    allv+=list(oi.values())+ (list(stta) if stta is not None else [])
ymin=np.floor((np.percentile(allv,0.5)-0.05)*20)/20; ymax=np.ceil((np.percentile(allv,99.5)+0.05)*20)/20
fig,axes=plt.subplots(NR,NC,figsize=(15.5,20),sharex=True,sharey=True); axes=axes.ravel()   # geometry matched to Fig S2
for k,l in enumerate(locs):
    ax=axes[k]; start,obs_y,oi,xf,tta,stta,sttap,prim=D[l]; fitmask=xf<=2019
    ax.axvspan(2020,2023,color="#f0e6d8",zorder=0); ax.axhline(0,color="#bbb",lw=.6); ax.axvline(2019,color="#bbb",lw=.6,ls=":")
    if tta is not None: ax.plot(xf[xf>=2019],tta[xf>=2019],color=CTT,lw=1.2,ls=(0,(3,2)),zorder=3)       # TTa projection
    if stta is not None:
        ax.plot(xf[fitmask],stta[fitmask],color=CST,lw=1.3,zorder=2)                                     # STTa fit
        ax.plot(xf[xf>=2019],stta[xf>=2019],color=CST,lw=1.3,ls=(0,(3,2)),zorder=2)                      # STTa projection
    if sttap is not None: ax.plot(sttap[0],sttap[1],color=CSP,lw=1.9,zorder=4)                           # STTa+ anchored curve
    if prim is not None:
        px,pv=prim
        ax.plot(px[px<=2019],pv[px<=2019],color=CBW,lw=1.5,zorder=6)
        ax.plot(px[px>=2019],pv[px>=2019],color=CBW,lw=1.7,ls=(0,(3,2)),zorder=6)
    ax.plot(obs_y,[oi[y] for y in obs_y],color=COBS,lw=1.5,zorder=5)                                     # observed
    ax.set_title(name.get(l,l),fontsize=8,fontweight="bold",pad=1.5); ax.tick_params(labelsize=6.3)   # matched to Fig S2
    ax.set_xlim(2003,2025); ax.set_ylim(ymin,ymax); ax.grid(alpha=.15); ax.set_axisbelow(True)
    ax.set_xticks([2005,2010,2015,2020,2025])                                                            # matched to Fig S2
    for s in ["top","right"]: ax.spines[s].set_visible(False)
for j in range(len(locs),NR*NC): axes[j].axis("off")
hp=[plt.Line2D([0],[0],color=COBS,lw=2.4,label="Observed ln(rate)"),
    plt.Line2D([0],[0],color=CTT,lw=2.4,ls=(0,(3,2)),label="All-ages two-step trend, own"),
    plt.Line2D([0],[0],color=CST,lw=2.4,label="All-ages two-step trend, shrunk"),
    plt.Line2D([0],[0],color=CSP,lw=2.4,label="All-ages, slope-anchored"),
    plt.Line2D([0],[0],color=CBW,lw=2.4,label="Primary STTa: standardised aggregate of the 20 band baselines"),
    plt.Line2D([0],[0],color="#f0e6d8",lw=10,label="2020–2023 (pandemic)")]
# legend starts at the left edge of the FIRST empty cell and grows rightwards across the
# remaining ones, instead of being centred in the last cell and spilling past the grid
axes[len(locs)].legend(handles=hp,loc="center left",bbox_to_anchor=(0.0,0.5),
                       fontsize=9.5,frameon=False,borderaxespad=0.0,handlelength=2.4)   # 9.5 = Fig S2
fig.supylabel("ln( age-standardised death rate ) − fitted 2019 level  (common scale)",fontsize=12.5)
# identical to Fig S2 so the two montages are a matched pair: same panel size, same positions
fig.subplots_adjust(left=0.07,right=0.995,top=0.975,bottom=0.03,wspace=0.06,hspace=0.25)
FIG=os.path.join(OUTD,"fig_option1_logquad_montage_v2.png"); fig.savefig(FIG,dpi=140,bbox_inches="tight"); plt.close(fig)
print("wrote",FIG,f"({len(locs)} populations; common y-axis {ymin:.2f}..{ymax:.2f})")
