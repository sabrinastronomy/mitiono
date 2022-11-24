import numpy as np
import corner
import matplotlib.pyplot as plt
import os
autocorr_only = True
truths = [10.274034, 21.226579, 0.1e-7, 0.1e-7]

min_ra = 10.2698673
min_dec = 21.222412
max_ra = 10.2782007
max_dec = 21.230746

lim_list = [[min_ra,
             max_ra],
            [min_dec,
             max_dec],
            [-5e-7, 5e-7],
            [-5e-7, 5e-7]]

# lim_list = [[1.02740342e+01 - 0.05 / 60 / 60,
#              1.02740342e+01 + 0.05 / 60 / 60],
#             [2.1226579e+01 - 0.25 / 60 / 60,
#              2.1226579e+01 + 0.25 / 60 / 60],
#             [-.8e-7, 1.2e-7],
#             [-.8e-7, 1.2e-7]]
ax_labels=['R.A. (deg)', 'Dec (deg)', r'$\Delta DM_{\rm CA}$',
                          r'$\Delta DM_{\rm CT}$']
dir_main = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/Current Research/mitiono/fringe fitting/important_results/"
plot_dir = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/Current Research/mitiono/fringe fitting/trace_plots/"
for file in os.listdir(dir_main):
    if file.endswith("50000_5e-07.npy"):
        full_file = os.path.join(dir_main, file)
        print(os.path.join(dir_main, file))
    else:
        continue
    res = np.load(full_file, allow_pickle=True).item()

    if autocorr_only:
        fname = file.split(".")[0]
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        var = 0
        fig.suptitle(fname)
        for i in range(2):
            for j in range(2):
                num_samples = len(res.samples.T[var])
                axes[i][j].hlines(truths[var], 0, num_samples, color="k", ls="--", zorder=1000, alpha=1, label="truth")
                axes[i][j].scatter(np.arange(num_samples), res.samples.T[var], c="g", s=0.01)
                axes[i][j].set_xlabel("sample number")
                axes[i][j].set_ylabel(ax_labels[var])
                axes[i][j].ticklabel_format(useOffset=False)
                axes[i][j].set_ylim(lim_list[var])
                var += 1

        plt.tight_layout()
        fig.savefig(plot_dir + f"{fname}_trace.png")

        # for i in range(len(truths)):
        #     num_samples = len(res.samples.T[i])
        #     plt.hlines(truths[i], 0, num_samples, color="k", ls="--", zorder=1000, alpha=1, label="truth")
        #     plt.title("30 arcseconds width, sigma_DM = 5e-7")
        #     plt.scatter(np.arange(num_samples), res.samples.T[i], c="g", s= 0.01)
        #     plt.savefig(dir_main + f"trace_{i}.png")
        #     plt.xlabel("sample number")
        #     plt.ylabel(ax_labels[i])
        #     plt.close()
        continue


    # print(np.shape(res))
    # plt.scatter(res[0])
    # plt.show()
    #
    # sys.exit()

    # plot initial run (res1; left)
    # Loop over the diagonal
    # print(np.shape(res))
    fig, ax = plt.subplots()


    truths = [10.274034, 21.226579, 0.1e-7, 0.1e-7]

    from matplotlib.ticker import ScalarFormatter, FuncFormatter
    def set_lim(fig, lim_list, axis='x'):
        ndim = int(np.sqrt(len(fig.axes)))
        for ii in range(ndim):
            for jj in range(ndim):
                if axis == 'x':
                    fig.axes[ndim * ii + jj].set_xlim(lim_list[jj])
                    x_formatter = ScalarFormatter(useOffset=False)
                    fig.axes[ndim * ii + jj].xaxis.set_major_formatter(x_formatter)
                if axis == 'y':
                    y_formatter = ScalarFormatter(useOffset=False)
                    fig.axes[ndim * ii + jj].yaxis.set_major_formatter(y_formatter)
                    if ii != jj:
                        fig.axes[ndim * ii + jj].set_ylim(lim_list[ii])

    def set_labels(fig, ax_labels):
        axs = fig.axes
        ndim = int(np.sqrt(len(axs)))
        for ii in range(ndim):
            for jj in range(ndim):
                # ii is row, jj is column
                if ii != 0 and jj == 0:  # left row
                    axs[ndim * ii + jj].set_ylabel(ax_labels[ii])
                    if ii == 1:
                        plt.setp(axs[ndim * ii + jj].get_yticklabels()[0],
                                 visible=False)  # remove first tick labels

                    if ii >= ndim - 2:  # DM axes
                        ticks_y = FuncFormatter(lambda x, pos: '{0:g}'.format(x / 1e-7))
                        axs[ndim * ii + jj].yaxis.set_major_formatter(ticks_y)

                    # axs[ndim*ii+jj].yaxis.set_major_locator(MaxNLocator(4,prune='both'))
                if ii == ndim - 1:  # bottom row
                    if jj == 1:
                        plt.setp(axs[ndim * ii + jj].get_xticklabels()[0], visible=False)  # remove tick labels
                    axs[ndim * ii + jj].set_xlabel(ax_labels[jj])
                    if jj >= ndim - 2:  # DM axes
                        ticks_x = FuncFormatter(lambda x, pos: '{0:g}'.format(x / 1e-7))
                        axs[ndim * ii + jj].xaxis.set_major_formatter(ticks_x)
                    # axs[ndim*ii+jj].xaxis.set_major_locator(MaxNLocator(4,prune='both'))

    def set_ticks_invisible(fig):
        axs = fig.axes
        ndim = int(np.sqrt(len(axs)))
        for ii in range(ndim):
            for jj in range(ndim):
                # ii is row, jj is column
                if jj > 0:
                    axs[ndim * ii + jj].axes.yaxis.set_ticklabels([])
                if ii < ndim - 1:
                    axs[ndim * ii + jj].axes.xaxis.set_ticklabels([])

    def flip_direction_column(fig, col=0, axis='x'):
        axs = fig.axes
        ndim = int(np.sqrt(len(axs)))
        for ii in range(ndim):
            for jj in range(ndim):
                # ii is row, jj is column
                if jj == col:
                    if axis == 'y':
                        axs[ndim * ii + jj].invert_yaxis()
                    if axis == 'x':
                        axs[ndim * ii + jj].invert_xaxis()

    mcmc_fig= corner.corner(res.samples, bins = 60, hist_bin_factor = 6, truths = truths)

    set_lim(mcmc_fig, lim_list, axis='x')
    set_lim(mcmc_fig, lim_list, axis='y')

    set_ticks_invisible(mcmc_fig)
    set_labels(mcmc_fig,
               ax_labels=['R.A. (deg)', 'Dec (deg)', r'$\Delta DM_{\rm CA} \times 10^{-7}$' + r'pc cm$^{-3}$',
                          r'$\Delta DM_{\rm CT} \times 10^{-7}$' + r'pc cm$^{-3}$'])  # r'$\Delta DM_{\rm DA} \times 10^{-7}\text{pc cm}^{-3}$', r'$\Delta DM_{\rm DG} \times 10^{-7}\text{pc cm}^{-3}$'])
    # (r'$\Delta DM_{\rm DG} \times 10^{-7}$'+ r'pc cm$^{-3}$')
    flip_direction_column(mcmc_fig, col=0, axis='x')

    plt.savefig(f"{plot_dir}{file}.pdf")
    plt.savefig(f"{plot_dir}{file}.png", dpi=200)