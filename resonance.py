import numpy as np

def draw_resonance_lines(ax, tunex_range=[0,1], tuney_range=[0,1], synchro_line=False,
                         synchro_order=[1,], synchro_freq=0.01,
                         max_order=8, line_width=0.4, coupling=True, working_point=None):
    n_dots = 2
    minimum_alpha=0.2
    maximum_alpha=0.6
    ax.set_aspect('equal')
    ax.set_xlabel('Horizontal Tune')
    ax.set_ylabel('Vertical Tune')
    ax.set_xlim(tunex_range)
    ax.set_ylim(tuney_range)

    if working_point is not None:
        ax.plot(working_point[0], working_point[1], marker='.')
    for order in range(1,max_order+1):
        xdots=np.linspace(tunex_range[0], tunex_range[1],n_dots)
        for l in range(1, order):
            ydots=np.zeros_like(xdots)+1.0*l/order
            ax.plot(xdots,ydots,'k--', alpha=minimum_alpha+maximum_alpha*(max_order-order)/max_order, linewidth=line_width)
            ax.plot(ydots,xdots,'k--', alpha=minimum_alpha+maximum_alpha*(max_order-order)/max_order, linewidth=line_width)
    if coupling==False:
        return
    for xorder in range(1,max_order):
        for yorder in range(1,max_order-xorder+1):
            order=xorder+yorder
            for l in range(-order+1, order):
                ydots=np.linspace(tuney_range[0], tuney_range[1],n_dots)
                xdots=(l-yorder*ydots)/xorder
                ax.plot(xdots,ydots,'k--', alpha=minimum_alpha+maximum_alpha*(max_order-abs(order))/max_order, linewidth=line_width*0.8)
                xdots=(l+yorder*ydots)/xorder
                ax.plot(xdots,ydots,'k--', alpha=minimum_alpha+maximum_alpha*(max_order-abs(order))/max_order, linewidth=line_width*0.8)
    if synchro_line:
        for sorder in synchro_order:
            for xorder in range(1, max_order-sorder):
                for yorder in range(1, max_order - sorder - xorder + 1):
                    order = xorder + yorder
                    for l in range(-order + 1, order):
                        ydots = np.linspace(tuney_range[0], tuney_range[1], n_dots)
                        xdots = (l - sorder*synchro_freq - yorder * ydots) / xorder
                        ax.plot(xdots, ydots, 'k--',
                                alpha=minimum_alpha + maximum_alpha * (max_order - abs(order)) / max_order,
                                linewidth=line_width * 0.8)
                        xdots = (l - sorder*synchro_freq + yorder * ydots) / xorder
                        ax.plot(xdots, ydots, 'k--',
                                alpha=minimum_alpha + maximum_alpha * (max_order - abs(order)) / max_order,
                                linewidth=line_width * 0.8)
                        xdots = (l + sorder * synchro_freq - yorder * ydots) / xorder
                        ax.plot(xdots, ydots, 'k--',
                                alpha=minimum_alpha + maximum_alpha * (max_order - abs(order)) / max_order,
                                linewidth=line_width * 0.8)
                        xdots = (l + sorder * synchro_freq + yorder * ydots) / xorder
                        ax.plot(xdots, ydots, 'k--',
                                alpha=minimum_alpha + maximum_alpha * (max_order - abs(order)) / max_order,
                                linewidth=line_width * 0.8)

    return
