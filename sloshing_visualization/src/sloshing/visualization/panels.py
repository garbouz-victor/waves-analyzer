"""Fixed-scale Matplotlib artists; updates change data, never normalization."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import numpy as np

from .scope import CONTACT_WARNING, CORNER_WARNING, FIXED_DOMAIN, magnify_positions, phase_label

plt.rcParams.update({"font.size": 11, "axes.titlesize": 13, "axes.labelsize": 11,
                     "figure.facecolor": "white", "savefig.facecolor": "white"})
LAYER_COLORS = {-0.1: "#f5b52e", -0.3: "#ea6f42", -1.: "#32bcad", -3.: "#a98ced"}


def display_cell_edges(nodes):
    """Half-cells at domain endpoints: no heatmap pixels outside z<=0.

    This limits the drawing geometry, not scalar values, eta or vorticity.
    """
    nodes=np.asarray(nodes)
    return np.r_[nodes[0],(nodes[:-1]+nodes[1:])/2,nodes[-1]]


def contact_mask(ax, scope, warning=True):
    artists = []
    a, b, c = (scope[k] for k in ("domain_abs_x_max_m", "bulk_abs_x_max_m", "transition_abs_x_max_m"))
    for side in (-1, 1):
        lo, hi = sorted((side*b, side*c))
        patch=Rectangle((lo,0),hi-lo,1,transform=ax.get_xaxis_transform(),color="#bd9f54",alpha=.09,zorder=5)
        ax.add_patch(patch);artists.append(patch)
        lo, hi = sorted((side*c, side*a))
        patch=Rectangle((lo,0),hi-lo,1,transform=ax.get_xaxis_transform(),facecolor=".7",edgecolor=".4",alpha=.6,hatch="////",lw=.6,zorder=8)
        ax.add_patch(patch);artists.append(patch)
    if warning:
        ax.text(.5, -.14, CONTACT_WARNING, ha="center", transform=ax.transAxes, fontsize=10, color=".3")
    return artists


def contact_corner_mask(ax,scope,warning=True):
    patches=[]
    low,high=scope['contact_corner_z_min_m'],scope['contact_corner_z_max_m']
    for lo,hi in ((-1.,-scope['transition_abs_x_max_m']),(scope['transition_abs_x_max_m'],1.)):
        p=Rectangle((lo,low),hi-lo,high-low,facecolor='.7',edgecolor='.4',alpha=.6,hatch='////',lw=.6,zorder=8)
        ax.add_patch(p);patches.append(p)
    if warning:ax.text(.5,-.14,CORNER_WARNING,ha='center',transform=ax.transAxes,fontsize=10,color='.3')
    return patches


def solid_walls(ax,zmin,bottom=False,sides=(-1,1)):
    segments=[([side,side],[zmin,0]) for side in sides]
    if bottom:segments.append(([-1,1],[-10,-10]))
    lines=[]
    for x,z in segments:
        ax.plot(x,z,color='white',lw=4.5,zorder=29,clip_on=False)
        line,=ax.plot(x,z,color='black',lw=2.3,zorder=30,clip_on=False,label='_nolegend_')
        lines.append(line)
    return lines


class SurfacePanel:
    def __init__(self, ax, h, meta, title="TRUE-SCALE FREE SURFACE"):
        self.ax, self.h = ax, h
        self.x = h["eta_x"][:]
        self.line, = ax.plot(self.x, h["eta"][0]*1000, color="#005e8a", lw=2)
        ax.plot(self.x, h["eta"][0]*1000, "--", color=".6", lw=1, label="initial")
        ax.axhline(0, color=".5", lw=.6)
        ax.scatter(self.x[[0,-1]], h["eta"][0][[0,-1]]*1000, c="black", s=22, zorder=10)
        contact_mask(ax, meta["scope"], warning=False)
        ax.set(xlim=(-1,1), ylim=meta["surface_ylim_mm"], xlabel="x (m)", ylabel="eta (mm)", title=title)
        ax.text(.98, .92, "pinned endpoints", ha="right", transform=ax.transAxes, fontsize=9)
        ax.text(.5,.04,"hatched contacts: unresolved",ha="center",transform=ax.transAxes,fontsize=8,color=".35")
        ax.grid(alpha=.2)

    def update(self, i):
        self.line.set_ydata(self.h["eta"][i]*1000)


class FlowPanel:
    def __init__(self, ax, h, meta, field="speed", tracers=True, full=False, path_key="linear_paths"):
        self.ax, self.h, self.meta = ax, h, meta
        self.field, self.full, self.path_key = field, full, path_key
        self.x = h["full_x" if full else "x"][:]
        self.z = h["full_z" if full else "z"][:]
        self.prefix = "full_" if full else ""
        limits = meta["velocity_color_range_m_per_s"] if field == "speed" else meta["vorticity_color_range_per_s"]
        self.image = ax.pcolormesh(display_cell_edges(self.x), display_cell_edges(self.z), np.zeros((len(self.z),len(self.x))), shading="flat",
                                  cmap="viridis" if field == "speed" else "RdBu_r", norm=Normalize(*limits), rasterized=True)
        label = "speed (m/s)" if field == "speed" else "omega = dw/dx - du/dz (1/s)"
        cb = ax.figure.colorbar(self.image, ax=ax, pad=.025, fraction=.034, extend="max" if field == "speed" else "both")
        cb.set_label(label)
        self.surface, = ax.plot(h["eta_x"][:], h["eta"][0], color="#146ea0", lw=1.5, zorder=10)
        ax.scatter([-1,1], h["eta"][0][[0,-1]], c="black", s=25, zorder=11)
        ax.axhline(0, color=".5", ls=":", lw=.8)
        self.masks = contact_corner_mask(ax, meta["scope"])
        ax.set(xlim=(-1,1), ylim=(-10 if full else -1.5,.04), xlabel="x (m)", ylabel="z (m)")
        ax.set_aspect("equal", adjustable="box")
        self.walls=solid_walls(ax,-10 if full else -1.5,bottom=full)
        self.wall_label=ax.annotate('no-slip: u = w = 0',xy=(-1,-1.08),xytext=(-.76,-1.10),
            fontsize=10,zorder=31,bbox={'facecolor':'white','alpha':.9,'edgecolor':'none','pad':3},
            arrowprops={'arrowstyle':'-','color':'black','lw':1.2})
        ax.text(.5,1.01,"surface x1; reference domain unchanged", ha="center", transform=ax.transAxes, fontsize=10)
        self.q = None
        if not full:
            self.q = ax.quiver(h["arrow_x"][:], h["arrow_z"][:], np.zeros_like(h["arrow_x"][:]), np.zeros_like(h["arrow_z"][:]),
                               color="white" if field == "speed" else "#232323", angles="xy", scale_units="xy",
                               scale=1/meta["fixed_arrow_seconds"], width=.0024, headwidth=3.4, pivot='tail', zorder=4)
            # Suppress Matplotlib's minimum-length dot at exact zero velocity.
            self.q.minlength = 0
            ax.quiverkey(self.q, .70, .045, meta["arrow_key_speed_m_per_s"], "0.5 mm/s (fixed arrows)",
                         labelpos="N", coordinates="axes", color="white" if field == "speed" else "black",
                         labelcolor="white" if field == "speed" else "black", fontproperties={"size":9})
        self.points = self.tails = None
        seeds = h["seeds"][:]
        self.selected = np.flatnonzero(seeds[:,1] >= -1.)
        if tracers:
            cols = [LAYER_COLORS[float(z)] for z in seeds[self.selected,1]]
            ax.scatter(seeds[self.selected,0], seeds[self.selected,1], c=cols, marker="+", s=22, alpha=.6, zorder=6)
            self.points = ax.scatter(seeds[self.selected,0], seeds[self.selected,1], c=cols, s=29, edgecolors="black", linewidths=.4, zorder=7)
            self.tails = LineCollection([], colors=cols, linewidths=1.3, alpha=.7, zorder=6)
            ax.add_collection(self.tails)
            for z, color in list(LAYER_COLORS.items())[:3]:
                ax.plot([],[],"o",color=color,label=f"z0={z:g} m",ms=5)
            ax.legend(loc="lower left", bbox_to_anchor=(.015,.115), ncol=3, framealpha=.85, fontsize=9)
            ax.text(.5,-.21, f"Linear tracer displacement x{meta['particle_displacement_magnification']:g} (same x/z); + = reference point",
                    ha="center", transform=ax.transAxes, fontsize=10)

    def update(self, i):
        if self.field == "speed":
            values = np.hypot(self.h[self.prefix+"u"][i], self.h[self.prefix+"w"][i])
        else:
            values = self.h[self.prefix+"omega"][i]
        self.image.set_array(values.ravel())
        self.surface.set_ydata(self.h["eta"][i])
        if self.q is not None:
            self.q.set_UVC(self.h["arrow_u"][i], self.h["arrow_w"][i])
        if self.points is not None:
            seeds = self.h["seeds"][:]
            p = magnify_positions(seeds, self.h[self.path_key][max(0,i-100):i+1], self.meta["particle_displacement_magnification"])
            self.points.set_offsets(p[-1,self.selected])
            self.tails.set_segments([p[:,j] for j in self.selected])


class MainFigure:
    def __init__(self, h, meta, clean=False, omega=False):
        self.h, self.meta = h, meta
        self.fig = plt.figure(figsize=(19.2,10.8), dpi=100)
        self.title = self.fig.text(.045,.965,"",fontsize=20,weight="bold")
        self.phase = self.fig.text(.045,.922,"",fontsize=13)
        self.fig.text(.045,.018, FIXED_DOMAIN, fontsize=11, color=".3")
        ax = self.fig.add_axes([.055,.20,.555,.66])
        self.flow = FlowPanel(ax,h,meta,field="omega" if omega else "speed",tracers=not omega)
        self.surface = SurfacePanel(self.fig.add_axes([.69,.78,.275,.13]),h,meta)
        self.energy_ax = self.fig.add_axes([.69,.57,.275,.13])
        self.times = h["times"][:]
        for key,col,label in (("kinetic_energy","#bf6532","K"),("potential_energy","#146ea0","P"),("total_energy","#202020","E")):
            self.energy_ax.plot(self.times,h[key][:]*1e7,color=col,label=label,lw=1.5)
        self.energy_ax.set(xlim=(0,5),ylim=(0,1.06*h["total_energy"][0]*1e7),xlabel="physical t (s)",ylabel="energy x10^-7 m^4/s^2")
        self.energy_ax.legend(ncol=3,fontsize=10,loc="upper right")
        self.energy_ax.grid(alpha=.2)
        self.energy_mark = self.energy_ax.axvline(0,color="#b5303c",lw=1.3)
        self.modal_ax = self.fig.add_axes([.69,.405,.275,.075])
        self.modal_ax.plot(self.times,h["modal_eta"][:]*1000,color="#584585",lw=1.5)
        self.modal_ax.axhline(0,color=".5",lw=.6)
        for t in meta["modal"]["zero_crossings_s"]:
            self.modal_ax.axvline(t,color=".75",lw=.6,ls=":")
        self.modal_ax.set(xlim=(0,5),xlabel="physical t (s)",ylabel="modal eta (mm)", title=f"Measured period ~{meta['modal']['estimated_period_s']:.3f} s")
        self.modal_mark = self.modal_ax.axvline(0,color="#b5303c",lw=1.3)
        # Dedicated full-depth PROFILE, never hiding the wall in the flow view.
        self.depth_ax = self.fig.add_axes([.69,.05,.275,.075])
        self.depth, = self.depth_ax.plot(h["depths"][:],np.ones(len(h["depths"]))*np.nan,color="#b74d32",lw=1.5)
        self.depth_ax.set(yscale="log", ylim=(meta["depth_log_display_min"],meta["depth_q_max"]*1.2),xlim=(-10,0))
        self.depth_ax.set_title("Full-depth penetration: Q(z)",fontsize=12)
        self.depth_ax.set_xlabel("z (m)",fontsize=10)
        self.depth_ax.set_ylabel("Q, log scale\nm^(3/2)/s",fontsize=9)
        self.depth_ax.tick_params(labelsize=9)
        self.depth_ax.grid(alpha=.2)
        self.depth_ax.text(.02,.76,"zero off log axis",transform=self.depth_ax.transAxes,fontsize=8)
        self.clean, self.omega = clean, omega
        self.wall_profile=None
        self.wall_text=None
        if not clean:
            from .no_slip_panels import WallProfilePanel
            self.wall_profile=WallProfilePanel(self.fig.add_axes([.69,.205,.275,.115]),h,meta,compact=True)
            self.wall_text=self.fig.text(.055,.882,'',fontsize=11)
        if clean:
            self.energy_ax.set_visible(False)
            self.modal_ax.set_visible(False)
            self.depth_ax.set_position([.69,.28,.275,.30])
        if omega:
            self.fig.text(.055,.087,"Color range based on |x|<=0.98; pointwise corner extrema excluded from scale.",fontsize=12)
            self.fig.text(.055,.055,"Integral wall-vorticity structure is coherent (~5.7% field difference); corner extrema unresolved.",fontsize=11)

    def update(self,i):
        t = float(self.h["times"][i])
        self.flow.update(i)
        self.surface.update(i)
        self.energy_mark.set_xdata([t,t])
        self.modal_mark.set_xdata([t,t])
        if self.wall_profile is not None:
            self.wall_profile.update(i)
            l,r=(float(self.h['wall_'+s+'_speed_max'][i]) for s in ('left','right'))
            self.wall_text.set_text(f'FEM wall velocity (actual P2): left {l:.2g}, right {r:.2g} m/s — no-slip')
        q = self.h["depth_q"][i]
        self.depth.set_ydata(np.where(q>0,q,np.nan))
        c=self.meta["parameters"]
        self.title.set_text(f"{'Vorticity' if self.omega else 'Bulk viscous sloshing'}  |  t = {t:.3f} s  |  nu = {c['nu']:g} m²/s  |  alpha = {c['alpha_deg']:g}°")
        self.phase.set_text("Conditionally qualified bulk model; hatched contacts excluded" if self.clean else phase_label(t,self.meta["modal"]))
        return self.fig


class SurfaceFigure:
    def __init__(self,h,meta):
        self.h,self.meta=h,meta
        self.fig=plt.figure(figsize=(19.2,10.8),dpi=100)
        self.panel=SurfacePanel(self.fig.add_axes([.08,.25,.85,.60]),h,meta,"TRUE-SCALE eta in millimetres — no amplitude magnification")
        self.time=self.fig.text(.08,.91,"",fontsize=20)
        self.fig.text(.08,.05,FIXED_DOMAIN,fontsize=12)
        self.fig.text(.5,.12,CONTACT_WARNING,ha="center",fontsize=12)
    def update(self,i):
        self.panel.update(i)
        self.time.set_text(f"t = {self.h['times'][i]:.3f} s  |  "+phase_label(float(self.h['times'][i]),self.meta["modal"]))
        return self.fig


class TracerFigure:
    def __init__(self,h,meta):
        self.h,self.meta=h,meta
        self.fig=plt.figure(figsize=(19.2,10.8),dpi=100)
        self.time=self.fig.text(.07,.94,"",fontsize=19)
        self.difference=self.fig.text(.07,.895,"",fontsize=12)
        self.lines=[]
        self.seeds=h["seeds"][:]
        colors=[LAYER_COLORS[float(z)] for z in self.seeds[:,1]]
        for pos,title,key in (([.06,.25,.40,.60],"Linear first-order displacement","linear_paths"),([.55,.25,.40,.60],"Advected path in linear Eulerian field","advected_paths")):
            ax=self.fig.add_axes(pos)
            ax.set(xlim=(-1,1),ylim=(-3.2,.05),xlabel="x (m)",ylabel="z (m)",title=title)
            ax.set_aspect("equal",adjustable="box")
            contact_corner_mask(ax,meta["scope"],warning=False)
            solid_walls(ax,-3.2)
            ax.text(.5,.03,'Solid walls: u = w = 0',transform=ax.transAxes,ha='center',fontsize=9)
            ax.scatter(self.seeds[:,0],self.seeds[:,1],c=colors,marker="+",s=28)
            pts=ax.scatter(self.seeds[:,0],self.seeds[:,1],c=colors,s=38,edgecolors="black",lw=.5)
            trails=LineCollection([],colors=colors,lw=1.2)
            ax.add_collection(trails)
            surf,=ax.plot(h["eta_x"][:],h["eta"][0],color="#146ea0")
            ax.scatter([-1,1],h["eta"][0][[0,-1]],c="black",s=20)
            self.lines.append((key,pts,trails,surf))
        self.fig.text(.5,.18,CORNER_WARNING,ha="center",fontsize=11)
        self.fig.text(.5,.13,f"Tracer displacement x{meta['particle_displacement_magnification']:g} in BOTH panels / BOTH coordinates; surface x1",ha="center",fontsize=14)
        self.fig.text(.5,.073,"Difference is higher-order in perturbation amplitude. Do not interpret accumulated nonlinear drift physically.",ha="center",fontsize=13)
        self.fig.text(.5,.03,"Color = initial depth; + = initial reference point. These are displacement diagrams, not moving-mesh CFD.",ha="center",fontsize=12)
    def update(self,i):
        self.time.set_text(f"Tracer model comparison  |  t = {self.h['times'][i]:.3f} s")
        difference=float(np.linalg.norm(self.h['linear_paths'][i]-self.h['advected_paths'][i],axis=-1).max())
        self.difference.set_text(f"Maximum true linear/pathline separation at this time: {difference*1e6:.4f} micrometres (not magnified)")
        for key,pts,trails,surf in self.lines:
            p=magnify_positions(self.seeds,self.h[key][max(0,i-100):i+1],self.meta["particle_displacement_magnification"])
            pts.set_offsets(p[-1]);trails.set_segments([p[:,j] for j in range(len(self.seeds))])
            surf.set_ydata(self.h["eta"][i])
        return self.fig
