import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rc('text', usetex=True)


def plot(x,y, z=None,**kwargs):
    """_summary_

    Args:
        x (_type_): _description_
        y (_type_): _description_
        z (_type_, optional): _description_. Defaults to None.
        ax (_type_, optional): _description_. Defaults to None.
    """
    fig = plt.figure();
    ax = fig.add_subplot();
    ax.margins(x = 0.00) # 5% padding in all directions
    ax.tick_params(axis='both', which='major', labelsize=16);

    ts =1
    key="textscale";
    if( key in kwargs  ):  
        ts = kwargs[key];
        
    label=""
    key="label";
    if( key in kwargs  ):  
        label = kwargs[key];

    labelfs =int(18*ts);

    key="xlabel";
    if( key in kwargs  ):  
        xlabel = kwargs[key];
        ax.set_xlabel( r"$ "+xlabel+" $", fontsize=labelfs );

    key="ylabel";
    if( key in kwargs  ):  
        ylabel = kwargs[key];
        ax.set_ylabel( r"$ "+ylabel+" $", fontsize=labelfs);
        
#    if z is not None:
#        vmin = np.min(z);
#        vmax = np.max(z);        
#        im   = _ax.scatter(x,y,s=1, c=z,cmap=None,vmin=vmin, vmax=vmax);
#        fig.colorbar(im, ax=_ax);
    
    ax.plot(x,y, label=label);


def add_legends( ax, legends ):
    
    lines = [ line for line in ax.lines ]
    for i,legend in enumerate(legends):
        lines[i].set_label(r"$"+legend+" $")
        ax.legend(fontsize=int(18), frameon=False,labelspacing=0.2, handlelength=1,borderpad=0.2,handletextpad=0.2);
        
        
        
def plot4fig( Xs, Ys, **kwargs ):
    #Since this will plot 4 figure, it is assume Ys has len 4
    num_figs = 4;
    assert len(Ys) == num_figs, "The dataset submited through Ys should have length 4"
    Ys = np.array(Ys);
    
    
    key = "shareX";
    sharex = kwargs[key] if key in kwargs else False;
    if sharex:
        skey = "Xlabels";
        if skey in kwargs:
            labels = kwargs[skey];
            labels= np.broadcast_to( labels, (num_figs, len(labels)) );
        Xs    = np.broadcast_to( Xs    , (num_figs, len(Xs)) )

    key = "fig";
    fig, axs = None, None;
    if key in kwargs:
        fig = kwargs[key];
        axs = fig.axes;
    else:
        fig, axs = plt.subplots( 2, 2, dpi=300 );
        axs = axs.flatten();

    key = "Ylabels";
    if key in kwargs:
        for (ax,label) in zip(axs,kwargs[key]):
            ax.set_ylabel(label);

    key = "Xlabels";
    if key in kwargs:
        for (ax,label) in zip(axs,kwargs[key]):
            ax.set_xlabel(label);

    for ax,X,Ylist in zip(axs,Xs,Ys):
        Ylist = [Ylist,] if len(Ylist.shape) == 1 else Ylist;
        for Y in Ylist:
            ax.plot(X,Y);
        ax.margins(x = 0.00) # 5% padding in all directions
    
    plt.tight_layout(pad=0.5, h_pad=None, w_pad=None);
    return fig;

def plot_ZwithArrows(XX,YY,ZZ,UU,VV, Zlabel="", mask = None,npoints=(3,3), ax=None):
    if mask is not None:
        UU,VV,ZZ = [ np.ma.array(XX, mask=mask) for XX in (UU,VV,ZZ) ];

    # Interpolate to regularly-spaced quad grid.
    x,y = XX[:,0].flatten(),YY[0].flatten();
    fzi,fui,fvi =[ RegularGridInterpolator((x,y),  V) for V in [ZZ,UU,VV] ];

    #numdivisions
    kgrid= np.meshgrid(*[ np.linspace(x.min(),x.max(),n) for n,x in zip(npoints,[XX,YY]) ],indexing='ij');
    xi,yi= kgrid;
    kgrid= np.transpose(kgrid).swapaxes(0,1);
    
    zi,ui,vi =[ fx(kgrid) for fx in (fzi,fui,fvi)];

    if mask is not None:
        UU,VV,ZZ = [ np.ma.array(XX, mask=mask(ZZ)) for XX in (UU,VV,ZZ) ];
        ui,vi,zi = [ np.ma.array(xi, mask=mask(zi)) for xi in (ui,vi,zi) ];

    plotcb = False;
    if ax is  None:
        ax = plt.gca();
        plotcb = True;

    cos = ax.contour (XX,YY, ZZ,levels=[0], cmap="bone")
    cs = ax.contourf(XX,YY, ZZ,levels=100, cmap="bone")
    Qv = ax.quiver(xi,yi,ui,vi, color="C1",  scale=10.0, width=0.021/2.0);
    for c in cs.collections:
        c.set_edgecolor("face")
    for c in cos.collections:
        c.set_edgecolor("white")
    if plotcb:
        cbar = plt.gcf().colorbar(cs, shrink=0.9,label=Zlabel);
    return 0;

    def band_spin_texture(self, kpoints, bidx, uop = None, vop = None, zop=None, zlims = None, ax=None ):
        #Compute first the mean values using the model

        Ops= [uop, vop, zop ] ;

        def get_band_data( x, bidx ):
            y = None;
            if x.shape[0] == kpoints.shape[0]:
                y = x;
            else:
                y = self.compute_dispersion(kpoints, proj_op = x )[:,1];
            y = np.array(list(map(list,y))).T;
            return y[bidx];

        #If is a compatible array of data, add it directly to Ops
        Ops= [ get_band_data(op,bidx) for op in Ops ];

        uop, vop, zop = Ops;
        KX,KY,KZ = np.transpose(kpoints);

        # Create triangulation.
        triang = mtri.Triangulation(KX, KY)

        def get_mask( x, y, z ):
            return mtri.CubicTriInterpolator(triang, z, kind='geom')(x, y) > 0.;

        def xypoints( npts ):
            return np.meshgrid(np.linspace(KX.min(), KX.max(), npts), np.linspace(KY.min(), KY.max(), npts));

        # Set up the figure
        if( ax is None ):
            fig, ax = plt.subplots(dpi=400);

        # Interpolate to regularly-spaced quad grid.
        xi, yi = xypoints(npts=100)
        fzi    = mtri.CubicTriInterpolator(triang, zop, kind='geom');
        zi     = np.ma.array( fzi(xi,yi) , mask=get_mask(xi, yi,zop) )

        cs = ax.contourf(xi, yi, zi,levels=100, cmap="bone")
#        fig.colorbar(cs, ax=ax, shrink=0.9 , label=r"$\rm Energy (meV)$")
        # This is the fix for the white lines between contour levels
        for c in cs.collections:
            c.set_edgecolor("face")


        xi, yi = xypoints(npts=20)
        fui= mtri.CubicTriInterpolator(triang, uop, kind='geom');
        fvi= mtri.CubicTriInterpolator(triang, vop, kind='geom');
        ui = np.ma.array(fui(xi, yi), mask=get_mask(xi, yi,zop));
        vi = np.ma.array(fvi(xi, yi), mask=get_mask(xi, yi,zop));
        Q = ax.quiver(xi, yi, ui, vi, color="C1",  scale=15, width=0.022/4);

        return ax;
    