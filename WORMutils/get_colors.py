from .Error_Handler import Error_Handler

def get_colors(colormap, ncol, alpha=1.0, hide_traceback=True):

    """
    Get a list of rgb values for a matplotlib colormap
    
    Parameters
    ----------
    colormap : str
        Name of the colormap to color the scatterpoints. Accepts "WORM",
        "colorblind", or matplotlib colormaps.
        See https://matplotlib.org/stable/tutorials/colors/colormaps.html
        The "colorblind" colormap is referenced from Wong, B. Points of view:
        Color blindness. Nat Methods 8, 441 (2011).
        https://doi.org/10.1038/nmeth.1618
    
    ncol : int
        Number of colors to return in the list.
    
    alpha : float, default 1.0
        An alpha value between 0.0 (transparent) and 1.0 (opaque).

    hide_traceback : bool, default True
        Hide traceback message when encountering errors handled by this function?
        When True, error messages handled by this class will be short and to
        the point.
    
    Returns
    -------
    colors : list
        A list of rgb color tuples
    """
    
    err_handler = Error_Handler(clean=hide_traceback)
    
    qualitative_cmaps = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
                         'Dark2', 'Set1', 'Set2', 'Set3',
                         'tab10', 'tab20', 'tab20b', 'tab20c']
    
    if colormap == "colorblind":
        # colors from Wong B. 2011, https://doi.org/10.1038/nmeth.1618
        colors = [(0, 0, 0, alpha), # black
                  (230/255, 159/255, 0, alpha), # orange
                  (86/255, 180/255, 233/255, alpha), # sky blue
                  (0, 158/255, 115/255, alpha), # bluish green
                  (240/255, 228/255, 66/255, alpha), # yellow
                  (0, 114/255, 178/255, alpha), # blue
                  (213/255, 94/255, 0, alpha), # vermillion
                  (204/255, 121/255, 167/255, alpha)] # reddish purple
        if ncol <= len(colors):
            return colors[:ncol]
        else:
            print("Switching from 'colorblind' colormap to 'viridis' because there are {} variables to plot.".format(ncol))
            colormap = "viridis"
    elif colormap == "WORM":
        colors = [(0, 0, 0, alpha), # black
                  (22/255, 153/255, 211/255, alpha), # blue
                  (232/255, 86/255, 66/255, alpha), # red
                  (245/255, 171/255, 80/255, alpha), # orange
                  (115/255, 108/255, 168/255, alpha), # purple
                  (151/255, 208/255, 119/255, alpha), # green
                  (47/255, 91/255, 124/255, alpha), # dark blue
                  (119/255, 119/255, 119/255, alpha)] # gray
        if ncol <= len(colors):
            return colors[:ncol]
        else:
            print("Switching from 'WORM' colormap to 'viridis' because there are {} variables to plot.".format(ncol))
            colormap = "viridis"
            
    if colormap in qualitative_cmaps:
        # handle qualitative (non-continuous) colormaps
        colors = [plt.cm.__getattribute__(colormap).colors[i] for i in range(ncol)]
        colors = [(c[0], c[1], c[2], alpha) for c in colors]
    else:
        # handle sequential (continuous) colormaps
        norm = matplotlib.colors.Normalize(vmin=0, vmax=ncol-1)
        try:
            cmap = cm.__getattribute__(colormap)
        except:
            valid_colormaps = [cmap for cmap in dir(cm) if "_" not in cmap and cmap not in ["LUTSIZE", "MutableMapping", "ScalarMappable", "functools", "datad", "revcmap"]]
            err_handler.raise_exception("'{}'".format(colormap)+" is not a recognized matplotlib colormap. "
                    "Try one of these: {}".format(valid_colormaps))
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors = [m.to_rgba(i) for i in range(ncol)]
        colors = [(c[0], c[1], c[2], alpha) for c in colors]
    
    return colors