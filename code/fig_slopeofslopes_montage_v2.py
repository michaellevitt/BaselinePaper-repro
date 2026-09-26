"""
Figure 1, v2 (2026-09-19) — 38 panels of the age-standardised slope-of-change, with the trend the
PRIMARY model now uses overlaid.

Age-specific trends became the primary analysis on 19 Sep 2026, so the three country-level lines
this figure used to label TTa, STTa and STTa+ are no longer the models' parameters: each of the
twenty age bands now carries its own two-step trend.  Those three lines are still worth showing,
because the heterogeneity of the age-standardised trend across populations is what motivates the
whole approach, but they are relabelled as what they are, all-ages trends.

A fourth line is added: the baseline-deaths-weighted mean of the twenty band trends actually used
by the primary STTa, read from option1_country_params.csv (columns g19_bandwmean, sos_bandwmean).
Where it departs from the teal all-ages line, promoting age-specific trends changed that
population's baseline.  Figure S2 shows the twenty bands behind it.

Points are the empirical 5-year centred moving slopes of the ESP-2013 age-standardised log death
rate.  Slope line g(y) = g19 + sos*(y-2019); a rising line means a decelerating decline.

Reads output/option1_country_params.csv.  Writes output/fig_slopeofslopes_montage_v2.png .
"""
import os, sys, csv, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); OUTD=os.path.join(ROOT,"output")
sys.path.insert(0, HERE)
from excess_anchor_window import load, cell_for
FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
      "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
ESP={"0":1000,"1-4":4000,"5-9":5500,"10-14":5500,"15-19":5500,"20-24":6000,"25-29":6000,
     "30-34":6500,"35-39":7000,"40-44":7000,"45-49":7000,"50-54":7000,"55-59":6500,"60-64":6000,
     "65-69":5500,"70-74":5000,"75-79":4000,"80-84":2500,"85-89":1500,"90+":1000}; ESPsum=sum(ESP.values())
RELIABLE={"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,"POL":2008,"SVK":2006}
X0,XF,XP=2003,2019,2025          # display start, fit end, projection end
COL="#2E75B6"; CBW="#e67e22"   # deaths-weighted mean of band trends (primary)
CTT="#c0392b"; CST="#159090"; CSP="#6a3d9a"   # points / TTa / STTa / STTa+
data,name=load()
def lnR(loc,y):
    if any(cell_for(data,loc,y,a,"T") is None for a in FINE): return None
    r=sum(ESP[a]*cell_for(data,loc,y,a,"T")[0]/cell_for(data,loc,y,a,"T")[1] for a in FINE)/ESPsum
    return np.log(r) if r>0 else None
def moving_slopes(loc):
    series={y:lnR(loc,y) for y in range(X0-2,XF+1)}; out={}
    for c in range(X0,XF+1):
        win=[(y,series[y]) for y in range(c-2,c+3) if X0-2<=y<=XF and series.get(y) is not None]
        if len(win)>=3:
            xs=np.array([w[0] for w in win],float); ys=np.array([w[1] for w in win]); xc=xs-xs.mean()
            out[c]=100.0*(xc*(ys-ys.mean())).sum()/(xc**2).sum()
    return out
par={r["country"]:r for r in csv.DictReader(open(os.path.join(OUTD,"option1_country_params.csv")))}
def anch(loc):         # STTa+ slope-anchored line: (g19a, sosa) in %/yr, %/yr2 from the shrunk ANCHORED slope
    ba=par[loc].get("beta_anch_shrunk"); ga=par[loc].get("gamma_anch_shrunk")
    if not ba or not ga: return None
    return 100.0*float(ba), 200.0*float(ga)
def ret_point(loc):    # observed return slope through {2019,2024,2025}, plotted at its centroid (~2022.7) = the STTa+ anchor datum
    yrs=[y for y in (2019,2024,2025) if lnR(loc,y) is not None]
    if len(yrs)<2: return None
    xx=np.array(yrs,float); yy=np.array([lnR(loc,y) for y in yrs]); xc=xx-xx.mean()
    return xx.mean(), 100.0*(xc*(yy-yy.mean())).sum()/(xc**2).sum()
order=sorted(par, key=lambda l:par[l]["name"]); assert len(order)==38, f"expected 38, got {len(order)}"
slopes={l:moving_slopes(l) for l in order}
retpts={l:ret_point(l) for l in order}
# ---- STTa* : STTa refit through the return-slope star (figure display, centre convention) ----
def own_anch_fit(loc):                                  # OLS through reliable pre-pandemic centred slopes + the star
    d=slopes[loc]; rsl=RELIABLE.get(loc,X0); rp=retpts[loc]
    pts=[(x,d[x]) for x in d if x>=rsl]
    if rp is None or len(pts)<2: return None
    xs=np.array([x for x,_ in pts]+[rp[0]],float); ys=np.array([v for _,v in pts]+[rp[1]],float)
    xb=xs.mean(); b=((xs-xb)*(ys-ys.mean())).sum()/((xs-xb)**2).sum(); a=ys.mean()-b*xb
    return a+b*2019.0, b                                # (g19* at 2019, sos*)
ownA={l:own_anch_fit(l) for l in order}
_ok=[l for l in order if ownA[l] is not None]; _pw=np.array([float(par[l]["pop2019"]) for l in _ok])
g19G=np.average([ownA[l][0] for l in _ok],weights=_pw); sosG=np.average([ownA[l][1] for l in _ok],weights=_pw)
def sstar(loc):                                         # shrink own anchored fit toward global by w
    o=ownA[loc]
    if o is None: return None
    w=float(par[loc]["w"]); return w*o[0]+(1-w)*g19G, w*o[1]+(1-w)*sosG
sstarD={l:sstar(l) for l in order}
# common y-axis from all points + line endpoints (fit and projection)
allv=[v for d in slopes.values() for v in d.values()]
for l in order:
    p=par[l]
    if p["g19_pct"]:
        g,s=float(p["g19_pct"]),float(p["sos_pct"]); allv+=[g+s*(X0-2019),g,g+s*(XP-2019)]
    gs=100*float(p["beta_shrunk"]); ss=200*float(p["gamma_shrunk"]); allv+=[gs+ss*(X0-2019),gs,gs+ss*(XP-2019)]
    av=anch(l)
    if av is not None: allv+=[av[0],av[0]+av[1]*(XP-2019)]
allv+=[retpts[l][1] for l in order if retpts[l] is not None]
for l in order:
    if sstarD[l] is not None: allv+=[sstarD[l][0]+sstarD[l][1]*(X0-2019), sstarD[l][0]+sstarD[l][1]*(XP-2019)]
ymin=np.floor(np.percentile(allv,2)-0.4); ymax=np.ceil(np.percentile(allv,98)+0.4)

NC,NR=5,8
fig,axes=plt.subplots(NR,NC,figsize=(15.5,20),sharex=True,sharey=True); axes=axes.ravel()
for k,loc in enumerate(order):
    ax=axes[k]; p=par[loc]; rs=RELIABLE.get(loc,X0)
    ax.axhline(0,color="black",lw=.7,zorder=1); ax.axvspan(2020,2023,color="#f4d03f",alpha=.10,zorder=0)
    ax.axvline(2019,color="#999",lw=.7,ls=":",zorder=1)
    d=slopes[loc]; xs=sorted(d)
    solid=[(x,d[x]) for x in xs if x>=rs]; dash=[(x,d[x]) for x in xs if x<=rs]
    if dash: ax.plot([x for x,_ in dash],[v for _,v in dash],color=COL,lw=1.1,ls="--",marker="o",ms=2.2,zorder=3)
    ax.plot([x for x,_ in solid],[v for _,v in solid],color=COL,lw=1.5,marker="o",ms=2.4,zorder=3)
    yf=np.array([max(rs,X0),XF],float); yp=np.array([XF,XP],float)
    if p["g19_pct"]:            # TTa (own)
        g,s=float(p["g19_pct"]),float(p["sos_pct"])
        ax.plot(yf,g+s*(yf-2019),color=CTT,lw=2.0,zorder=5); ax.plot(yp,g+s*(yp-2019),color=CTT,lw=1.5,ls=(0,(3,2)),zorder=5)
    gs=100*float(p["beta_shrunk"]); ss=200*float(p["gamma_shrunk"])     # STTa (shrunk)
    ax.plot(yf,gs+ss*(yf-2019),color=CST,lw=1.9,zorder=6); ax.plot(yp,gs+ss*(yp-2019),color=CST,lw=1.5,ls=(0,(3,2)),zorder=6)
    av=anch(loc)                                                        # STTa* = engine's slope-anchored shrunk trend (≡ paper STTa+, centre convention)
    if av is not None:
        xx=np.array([max(rs,X0),XP],float); ax.plot(xx,av[0]+av[1]*(xx-2019),color=CSP,lw=2.2,zorder=7)
    gbw,sbw=p.get("g19_bandwmean"),p.get("sos_bandwmean")               # PRIMARY: deaths-weighted mean of the 20 band trends
    if gbw not in (None,"") and sbw not in (None,""):
        gb,sb=float(gbw),float(sbw)
        ax.plot(yf,gb+sb*(yf-2019),color=CBW,lw=2.0,zorder=8)
        ax.plot(yp,gb+sb*(yp-2019),color=CBW,lw=1.6,ls=(0,(3,2)),zorder=8)
    rp=retpts[loc]                                                      # observed return-slope anchor point at ~2022.7
    if rp is not None: ax.plot(rp[0],rp[1],marker="*",ms=11,color=CSP,mec="black",mew=0.6,ls="none",zorder=9)
    ttl=name.get(loc,loc) if loc not in RELIABLE else f"{name.get(loc,loc)} (from {RELIABLE[loc]})"
    ax.set_title(ttl,fontsize=8,fontweight="bold",color=("#8b0000" if loc in RELIABLE else "black"))
    ax.set_ylim(ymin,ymax); ax.set_xlim(X0-0.5,XP+0.5); ax.set_xticks([2005,2010,2015,2020,2025])
    ax.tick_params(labelsize=6.3); ax.grid(alpha=.3); ax.set_axisbelow(True)
for k in range(len(order),len(axes)): axes[k].axis("off")
handles=[plt.Line2D([],[],color=COL,lw=1.6,marker="o",ms=4,label="Observed pre-pandemic 5-year slopes"),
         plt.Line2D([],[],color=CTT,lw=2.2,label="All-ages two-step trend, own"),
         plt.Line2D([],[],color=CST,lw=2.0,label="All-ages two-step trend, shrunk"),
         plt.Line2D([],[],color=CSP,lw=2.4,label="All-ages, slope-anchored"),
         plt.Line2D([],[],color=CBW,lw=2.2,label="Age-band STTa: deaths-weighted mean of the 20 band trends"),
         plt.Line2D([],[],color=CSP,marker="*",ms=13,ls="none",mec="black",mew=0.6,label="Observed slope for 2019, 2024 and 2025 data points")]
# legend starts at the left edge of the FIRST empty cell and grows rightwards across the
# remaining ones, instead of being centred in the last cell and spilling past the grid (match Fig S1)
axes[len(order)].legend(handles=handles,loc="center left",bbox_to_anchor=(0.0,0.5),
                        fontsize=9.5,frameon=False,borderaxespad=0.0,handlelength=2.4)
fig.supylabel("annual slope of the ESP-2013 age-standardised death rate  (%/yr)",fontsize=12.5)
fig.subplots_adjust(left=0.07,right=0.995,top=0.975,bottom=0.03,wspace=0.06,hspace=0.25)   # manual layout (match Fig S1); tight_layout was shrinking panels around the in-cell legend
FIG=os.path.join(OUTD,"fig_slopeofslopes_montage_v2.png"); fig.savefig(FIG,dpi=140,bbox_inches="tight"); plt.close(fig)
print(f"wrote {FIG}  (y-axis common {ymin:.0f}..{ymax:.0f} %/yr; 38 panels)")
