"""Helper functions to export simulation results.
"""
import json
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from matplotlib import cm
from matplotlib import rcParams

from . import data

rcParams.update({'figure.figsize': (8, 6)})  # figure size in inches


def _write_data_file(cfg, file_name, values, label, exp=False):
    out_path = cfg['OUT'] / file_name
    out_calc_text = " ".join(label) + "\n"
    for out_data in zip(*values):
        if exp:
            out_calc_text += ' '.join([f'{v:.5E}' for v in list(out_data)]) + "\n"
        else:
            out_calc_text += ' '.join([f'{v:.5F}' for v in list(out_data)]) + "\n"
    out_path.write_text(out_calc_text)


def _get_factors(cfg, multi=False, segments=False):
    if segments:
        if cfg['MULTI']['segment_export'] != 'mass':
            vol_factor = 1 / cfg['CATALYST']['volume']
            surf_factor = 1 / cfg['CATALYST']['surface_total']
            return vol_factor, surf_factor
    if multi:
        key = 'molar_weight_multi'
    else:
        key = 'molar_weight'
    vol_factor = (cfg['MOLECULE'][key] / data.Parameter.avogadro) / cfg['CATALYST']['volume']
    surf_factor = (cfg['MOLECULE'][key] / data.Parameter.avogadro) / cfg['CATALYST']['surface_total']
    return vol_factor, surf_factor


def fit_dark(cfg, result, result_raw):
    """Export function for fit: *Adsorption - Desorption* experiment in the dark.

    :param cfg: simulation config
    :type cfg: dict
    :param result: fit_t, fit_y, k_ads, k_des
    :type result: tuple
    :param result_raw: raw fit function result

    :result:
        * plot with fit and initial data points ``fit_dark.pdf``
        * space separated data file containing the fit ``fit_dark-values_calculated.txt``
        * space separated data file containing initial data points ``fit_dark-values_experiment.txt``
        * .json with the fitted parameters ``fit_dark.json``
        * .txt files with corresponding units (LaTeX formatted)
    """
    fit_t, fit_y, k_ads, k_des = result
    vol_factor, surf_factor = _get_factors(cfg)

    file_path = cfg['OUT'] / 'fit_dark.json'
    exp_dic = dict(
        k_ads=f"{k_ads:.3E} m/s",
        k_des=f"{k_des:.3E} 1/s",
        error=f"{result_raw.fun:.3E}",
        iterations=result_raw.nit
    )
    file_path.write_text(json.dumps(exp_dic, indent=4))

    base_y = []
    if len(fit_y) == 1:
        base_y.append(surf_factor * (cfg['DATA']['initial_concentration'][0] + cfg['DATA']['initial_concentration'][1] - fit_y[0]))
    else:
        for n, y in enumerate(fit_y):
            base_y.append(surf_factor*(cfg['DATA']['initial_concentration'][0][n] + cfg['DATA']['initial_concentration'][1][n] - y))

    surf_scale, surf_unit = data.Parameter.scale_ten(np.max(base_y))
    t_scale, t_unit = data.Parameter.scale_time(np.max(fit_t))
    label_scale, label_unit = data.Parameter.scale_ten(np.max(cfg['DATA']['initial_concentration'][1]) * vol_factor)

    cmap = plt.get_cmap('copper')
    f = plt.figure()
    ax1 = f.add_subplot(111)
    ax1.set_ylabel(f"surface concentration $C$ ({surf_unit}g/m$^2$)")
    ax1.set_xlabel(f"time $t$ ({t_unit})")
    out_calc = [list(fit_t * t_scale)]

    for n, y in enumerate(base_y):
        plt.plot(fit_t * t_scale, surf_scale * y, c='xkcd:darkblue')
        out_calc.append(list(surf_scale * y))

    if cfg['DATA']['initial_concentration'].ndim == 2:
        out_exp = [list(cfg['DATA']['time_series'][0][0] * t_scale)]
        for n in range(cfg['DATA']['initial_concentration'].shape[1]):
            l_value = cfg['DATA']['initial_concentration'][1][n] * vol_factor * label_scale
            ax1.plot(cfg['DATA']['time_series'][n][0] * t_scale,
                     surf_factor * surf_scale * (cfg['DATA']['initial_concentration'][0][n]+cfg['DATA']['initial_concentration'][1][n]-cfg['DATA']['time_series'][n][1]),
                     'o', c=cmap(n/(cfg['DATA']['initial_concentration'].shape[1]-1)),
                     label=f'{l_value:4G}'.strip() + f' {label_unit}g/m$^3$')
            out_exp.append(list(surf_factor * surf_scale * (cfg['DATA']['initial_concentration'][0][n]+cfg['DATA']['initial_concentration'][1][n]-cfg['DATA']['time_series'][n][1])))
    else:
        out_exp = [list(cfg['DATA']['time_series'][0] * t_scale)]
        l_value = cfg['DATA']['initial_concentration'][1] * vol_factor * label_scale
        ax1.plot(cfg['DATA']['time_series'][0] * t_scale,
                 surf_factor * surf_scale * (
                             cfg['DATA']['initial_concentration'][0] + cfg['DATA']['initial_concentration'][1] -
                             cfg['DATA']['time_series'][1]),
                 'o', c=cmap(0),  label=f'{l_value:4G}'.strip() + f' {label_unit}g/m$^3$')
        out_exp.append(list(surf_factor * surf_scale * (
                             cfg['DATA']['initial_concentration'][0] + cfg['DATA']['initial_concentration'][1] -
                             cfg['DATA']['time_series'][1])))

    ax1.set_ylim(0, ax1.get_ylim()[1])
    lg = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, ncol=1,
                    title='initial\nconcentration')
    lg.get_frame().set_alpha(0.8)
    lg.get_frame().set_linewidth(0.5)
    f.tight_layout()
    f.savefig(cfg['OUT']/'fit_dark.pdf')
    plt.close(f)

    out_path = cfg['OUT'] / 'fit_dark-unit.txt'
    out_path.write_text(f'{surf_unit}g/m$^2$')

    out_path = cfg['OUT'] / 'fit_dark-unit_time.txt'
    out_path.write_text(f'{t_unit}')

    _write_data_file(cfg, 'fit_dark-values_calculated.txt',
                     values=out_calc, label=['t'] + [f'fit{n:02}' for n in range(len(out_calc)-1)])

    _write_data_file(cfg, 'fit_dark-values_experiment.txt',
                     values=out_exp, label=['t'] + [f'exp{n:02}' for n in range(len(out_exp)-1)])


def fit_single(cfg, result, result_raw):
    """Export function for fit: *single species* degradation model.

    :param cfg: simulation config
    :type cfg: dict
    :param result: fit_t, fit_y, k_ads, k_des, k_reac
    :type result: tuple
    :param result_raw: raw fit function result

    :result:
        * plot with fit and initial data points ``fit_single.pdf``
        * space separated data file containing the fit ``fit_single-values_calculated.txt``
        * space separated data file containing initial data points ``fit_single-values_experiment.txt``
        * .json with the fitted parameters ``fit_single.json``
        * .txt files with corresponding units (LaTeX formatted)
    """
    cmap = plt.get_cmap('copper')
    fit_t, fit_y, k_ads, k_des, k_reac = result

    vol_factor, surf_factor = _get_factors(cfg)

    file_path = cfg['OUT'] / 'fit_single.json'
    exp_dic = dict(
        k_ads=f"{k_ads:.3E} m/s",
        k_des=f"{k_des:.3E} 1/s",
        k_reac=f"{k_reac:.3E} 1/s",
        error=f"{result_raw.fun:.3E}",
        iterations=result_raw.nit
    )
    file_path.write_text(json.dumps(exp_dic, indent=4))

    base_y = []
    if len(fit_y) == 1:
        base_y.append(
            vol_factor * fit_y[0])
    else:
        for n, y in enumerate(fit_y):
            base_y.append(vol_factor * (y))

    vol_scale, vol_unit = data.Parameter.scale_ten(np.max(base_y))
    t_scale, t_unit = data.Parameter.scale_time(np.max(fit_t))

    f = plt.figure()
    ax1 = f.add_subplot(111)
    ax1.set_ylabel(f"volume concentration $C$ ({vol_unit}g/m$^3$)")
    ax1.set_xlabel(f"time $t$ ({t_unit})")
    out_calc = [list(fit_t * t_scale)]

    for n, y in enumerate(base_y):
        plt.plot(fit_t * t_scale, vol_scale * y, c='xkcd:darkblue')
        out_calc.append(list(vol_scale * y))

    if cfg['DATA']['time_series'].ndim == 3:
        out_exp = [list(cfg['DATA']['time_series'][0][0] * t_scale)]
        for n in range(cfg['DATA']['time_series'].shape[0]):
            ax1.plot(cfg['DATA']['time_series'][n][0] * t_scale,
                     vol_factor * vol_scale * (cfg['DATA']['time_series'][n][1]),
                     'o', c=cmap(n / (cfg['DATA']['time_series'].shape[0] - 1)),
                     label=f'EXP {n:02}')
            out_exp.append(list(vol_factor * vol_scale * (
                    cfg['SYSTEM']['concentration_solution'] - cfg['DATA']['time_series'][n][1])))
    else:
        out_exp = [list(cfg['DATA']['time_series'][0] * t_scale)]
        ax1.plot(cfg['DATA']['time_series'][0] * t_scale,
                 vol_factor * vol_scale * cfg['DATA']['time_series'][1],
                 'o', c=cmap(0),
                 label='EXP 01')
        out_exp.append(list(vol_factor * vol_scale * (
                cfg['SYSTEM']['concentration_solution'] - cfg['DATA']['time_series'][1])))

    ax1.set_ylim(0, ax1.get_ylim()[1])
    lg = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, ncol=1)
    lg.get_frame().set_alpha(0.8)
    lg.get_frame().set_linewidth(0.5)
    f.tight_layout()
    f.savefig(cfg['OUT'] / 'fit_single.pdf')
    plt.close(f)

    out_path = cfg['OUT'] / 'fit_single-unit.txt'
    out_path.write_text(f'{vol_unit}g/m$^3$')

    out_path = cfg['OUT'] / 'fit_single-unit_time.txt'
    out_path.write_text(f'{t_unit}')

    _write_data_file(cfg, 'fit_single-values_calculated.txt',
                     values=out_calc, label=['t'] + [f'fit{n:02}' for n in range(len(out_calc)-1)])

    _write_data_file(cfg, 'fit_single-values_experiment.txt',
                     values=out_exp, label=['t'] + [f'exp{n:02}' for n in range(len(out_exp)-1)])


def fit_toc(cfg, sd_error, t, fit):
    """Export function for fit: *TOC*.

    :param cfg: simulation config
    :type cfg: dict
    :param sd_error: standard derivation error
    :type sd_error: float
    :param t: fit time series
    :type t: ndarray
    :param fit: fit toc series
    :type fit: ndarray

    :result:
        * plot with fit and initial data points | ``fit_toc.pdf``
        * space separated data file containing the fit | ``fit_toc-values_calculated.txt``
        * space separated data file containing initial data points | ``fit_toc-values_experiment.txt``
        * .json with the fitted parameters | ``fit_toc.json``
        * .txt files with corresponding units (LaTeX formatted)
    """
    exp_t = cfg['DATA']['time_series'][0]
    exp_data = cfg['DATA']['time_series'][1]
    fit = exp_data[0] * fit

    toc_factor = (12.0107 / data.Parameter.avogadro) / cfg['CATALYST']['volume']
    toc_scale, toc_unit = data.Parameter.scale_ten(exp_data[0]*toc_factor)
    t_scale, t_unit = data.Parameter.scale_time(np.max(t))

    f = plt.figure()
    ax1 = f.add_subplot(111)

    ax1.set_ylabel(f"TOC ({toc_unit}g/m$^3$)")
    ax1.set_xlabel(f"time $t$ ({t_unit})")

    ax1.plot(t_scale * t, fit*toc_scale*toc_factor, c='xkcd:blue', label='simulation')
    ax1.plot(t_scale * exp_t, exp_data*toc_scale*toc_factor, c='xkcd:orange', marker='o', linestyle='None',
             label='experiment')

    ax1.set_ylim(0, ax1.get_ylim()[1])
    ax1.set_xlim(0, ax1.get_xlim()[1])

    lg = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, ncol=1)
    lg.get_frame().set_alpha(0.8)
    lg.get_frame().set_linewidth(0.5)
    f.tight_layout()
    f.savefig(cfg['OUT'] / 'fit_toc.pdf')
    plt.close(f)

    file_path = cfg['OUT'] / 'fit_toc.json'
    exp_dic = dict(
        k_ads=f"{cfg['SYSTEM']['k_ads']:.3E} m/s",
        k_des=f"{cfg['SYSTEM']['k_des']:.3E} 1/s",
        k_reac=f"{cfg['SYSTEM']['k_reac']:.3E} 1/s",
        beta_0=f"{cfg['MULTI_WEAK']['beta_0']:.3E} 1/s",
        beta_1=f"{cfg['MULTI_WEAK']['beta_1']:.3E} 1/s",
        sd_error=f"{sd_error:.3E}"
    )
    file_path.write_text(json.dumps(exp_dic, indent=4))

    _write_data_file(cfg, 'fit_toc-values_calculated.txt',
                     values=[t*t_scale, fit*toc_scale*toc_factor], label=['t', 'fit'])

    _write_data_file(cfg, 'fit_toc-values_experiment.txt',
                     values=[exp_t*t_scale, exp_data*toc_scale*toc_factor], label=['t', 'exp'])

    out_path = cfg['OUT'] / 'fit_toc-unit.txt'
    out_path.write_text(f'{toc_unit}g/m$^3$')
    out_path = cfg['OUT'] / 'fit_toc-unit_time.txt'
    out_path.write_text(f'{t_unit}')


def single_species(cfg, t, N_surf, N_vol):
    """Export function single species model simulations.

        :param cfg: simulation config
        :type cfg: dict
        :param t: fit time series
        :type t: ndarray
        :param N_surf: number of molecules on the surface
        :type N_surf: ndarray
        :param N_vol: number of molecules in solution
        :type N_vol: ndarray

        :result:
            * plot with concentration development in solution and on the surface ``single_species.pdf``
            * space separated data file containing the concentrations in solution ``single_species-values_volume.txt``
            * space separated data file containing the concentrations on the surface ``single_species-values_surface.txt``
            * .txt files with corresponding units (LaTeX formatted)
        """

    vol_factor, surf_factor = _get_factors(cfg)

    vol_scale, vol_unit = data.Parameter.scale_ten(np.max(N_vol) * vol_factor)
    surf_scale, surf_unit = data.Parameter.scale_ten(np.max(N_surf) * surf_factor)
    t_scale, t_unit = data.Parameter.scale_time(np.max(t))

    f = plt.figure()
    ax1 = f.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.set_ylabel(f"volume concentration $C$ ({vol_unit}g/m$^3$)", color='xkcd:blue')
    ax2.set_ylabel(f"surface concentration $C$ ({surf_unit}g/m$^2$)", color='xkcd:orange')
    ax1.set_xlabel(f"time $t$ ({t_unit})")
    
    ax1.plot(t_scale * t, N_vol * vol_factor * vol_scale, c='xkcd:blue')
    ax2.plot(t_scale * t, N_surf * surf_factor * surf_scale, c='xkcd:orange')

    ax1.set_ylim(0, ax1.get_ylim()[1])
    ax2.set_ylim(0, ax2.get_ylim()[1])

    f.tight_layout()
    f.savefig(cfg['OUT'] / 'single_species.pdf')
    plt.close(f)

    out_path = cfg['OUT'] / 'single_species-unit_volume.txt'
    out_path.write_text(f'{vol_unit}g/m$^3$')

    out_path = cfg['OUT'] / 'single_species-unit_surface.txt'
    out_path.write_text(f'{surf_unit}g/m$^2$')

    out_path = cfg['OUT'] / 'single_species-unit_time.txt'
    out_path.write_text(f'{t_unit}')

    _write_data_file(cfg, 'single_species-values_volume.txt',
                     values=[t_scale * t, vol_scale * vol_factor * N_vol], label=['t', 'vol'])

    _write_data_file(cfg, 'single_species-values_surface.txt',
                     values=[t_scale * t, surf_scale * surf_factor * N_surf], label=['t', 'surf'])


def multi_species(cfg, t, N_surf, N_vol):
    """Export function for *incremental* and *fragmentation* multi species model simulations.

    :param cfg: simulation config
    :type cfg: dict
    :param t: fit time series
    :type t: ndarray
    :param N_surf: number of molecules on the surface
    :type N_surf: ndarray
    :param N_vol: number of molecules in solution
    :type N_vol: ndarray

    :result:
        * plot with concentration development in solution and TOC ``c_volume-toc.pdf``
        * plot of the segments in solution and on the surface ``volume_segments.pdf``, ``solution_segments.pdf``
        * space separated data file containing the concentrations in solution ``multi_species-values_volume.txt``
        * space separated data file containing the concentrations on the surface ``multi_species-values_surface.txt``
        * .txt files with corresponding units (LaTeX formatted)
    """

    c_count = len(N_surf[0])
    line_cm = cm.get_cmap('plasma_r', c_count - 2)
    norm = matplotlib.colors.Normalize(vmin=1.5, vmax=c_count - 0.5)
    vol_factor, surf_factor = _get_factors(cfg, multi=True)
    vol_seg_factor, surf_seg_factor = _get_factors(cfg, multi=True, segments=True)
    toc_factor = (12.0107 / data.Parameter.avogadro) / cfg['CATALYST']['volume']

    # Volume segments
    f = plt.figure()
    ax1 = f.add_subplot(111)
    vol_scale, vol_unit = data.Parameter.scale_ten(np.max(N_vol * vol_factor))
    t_scale, t_unit = data.Parameter.scale_time(np.max(t))
    if cfg['MULTI']['segment_export'] == 'mass':
        ax1.set_ylabel(f"volume concentration $C$ ({vol_unit}g/m$^3$)")
    else:
        ax1.set_ylabel(f"volume concentration $C$ (1/m$^3$)")
    ax1.set_xlabel(f"time $t$ ({t_unit})")
    ax1.set_yscale('log')
    for k in range(1, c_count-1):
        if cfg['MULTI']['segment_export'] == 'mass':
            ax1.plot(t_scale * t, N_vol[:, k] * vol_seg_factor[k] * vol_scale, c=line_cm(norm(k + 1)))
        else:
            ax1.plot(t_scale * t, N_vol[:, k] * vol_seg_factor, c=line_cm(norm(k + 1)))
    ax1.set_xlim(0, np.max(t_scale*t))
    calculated_ylim = ax1.get_ylim()
    ax1.set_ylim((np.max((calculated_ylim[0], 0.1 * vol_factor[-1] * vol_scale)), calculated_ylim[1]))
    cbar = f.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=line_cm), orientation='vertical', label='carbon count')
    cbar.ax.invert_yaxis()
    f.tight_layout()
    f.savefig(cfg['OUT'] / 'volume_segments.pdf')
    plt.close(f)

    # Surface segments
    f = plt.figure()
    ax1 = f.add_subplot(111)
    surf_scale, surf_unit = data.Parameter.scale_ten(np.max(N_surf * surf_factor))
    t_scale, t_unit = data.Parameter.scale_time(np.max(t))
    if cfg['MULTI']['segment_export'] == 'mass':
        ax1.set_ylabel(f"surface concentration $C$ ({surf_unit}g/m$^2$)")
    else:
        ax1.set_ylabel(f"surface concentration $C$ (1/m$^2$)")

    ax1.set_xlabel(f"time $t$ ({t_unit})")
    ax1.set_yscale('log')
    for k in range(1, c_count-1):
        if cfg['MULTI']['segment_export'] == 'mass':
            ax1.plot(t_scale * t, N_surf[:, k] * surf_seg_factor[k] * surf_scale, c=line_cm(norm(k + 1)))
        else:
            ax1.plot(t_scale * t, N_surf[:, k] * surf_seg_factor, c=line_cm(norm(k + 1)))
    ax1.set_xlim(0, np.max(t_scale*t))
    calculated_ylim = ax1.get_ylim()
    ax1.set_ylim((np.max((calculated_ylim[0], 0.1 * surf_factor[-1] * surf_scale)), calculated_ylim[1]))
    cbar = f.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=line_cm), orientation='vertical', label='carbon count')
    cbar.ax.invert_yaxis()
    f.tight_layout()
    f.savefig(cfg['OUT'] / 'surface_segments.pdf')
    plt.close(f)

    # Initial molecule and TOC
    f = plt.figure()
    ax1 = f.add_subplot(111)
    ax2 = ax1.twinx()
    toc = N_vol[:, 1:]
    if cfg['MULTI']['TOC_estimation'] == 'all':
        toc += N_surf[:, 1:]
    toc = np.sum(toc * (np.arange(c_count-1) + 2), axis=1)
    toc_scale, toc_unit = data.Parameter.scale_ten(toc[0] * toc_factor)
    vol_scale, vol_unit = data.Parameter.scale_ten(np.max(N_vol * vol_factor[-1]))
    t_scale, t_unit = data.Parameter.scale_time(np.max(t))
    ax1.set_ylabel(f"volume concentration $C$ ({vol_unit}g/m$^3$)", color='xkcd:blue')
    ax2.set_ylabel(f"TOC $C$ ({toc_unit}g/m$^3$)", color='xkcd:orange')
    ax1.set_xlabel(f"time $t$ ({t_unit})")
    ax1.plot(t_scale * t, N_vol[:, -1] * vol_factor[-1] * vol_scale, c='xkcd:blue')
    ax2.plot(t_scale * t, toc * toc_factor * toc_scale, c='xkcd:orange')
    ax1.set_ylim(0, ax1.get_ylim()[1])
    ax2.set_ylim(0, ax2.get_ylim()[1])

    f.tight_layout()
    f.savefig(cfg['OUT'] / 'c_volume-toc.pdf')
    plt.close(f)

    out_path = cfg['OUT'] / 'multi_species-unit_time.txt'
    out_path.write_text(f'{t_unit}')

    out_path = cfg['OUT'] / 'multi_species-unit_toc.txt'
    out_path.write_text(f'{toc_unit}g/m^3')

    _write_data_file(cfg, 'multi_species-values_toc.txt',
                     values=[t_scale * t, toc * toc_factor * toc_scale], label=['t', 'toc'])

    if cfg['MULTI']['segment_export'] == 'mass':
        out_path = cfg['OUT'] / 'multi_species-unit_volume.txt'
        out_path.write_text(f'{vol_unit}g/m$^3$')

        out_path = cfg['OUT'] / 'multi_species-unit_surface.txt'
        out_path.write_text(f'{surf_unit}g/m$^2$')

        _write_data_file(cfg, 'multi_species-values_volume.txt',
                         values=np.hstack((t_scale * t[:, np.newaxis], vol_scale * vol_factor * N_vol)).T,
                         label=['t'] + [f'v{n+1:02}' for n in range(c_count)])

        _write_data_file(cfg, 'multi_species-values_surface.txt',
                         values=np.hstack((t_scale * t[:, np.newaxis], surf_scale * surf_factor * N_surf)).T,
                         label=['t'] + [f's{n + 1:02}' for n in range(c_count)])

    else:
        out_path = cfg['OUT'] / 'multi_species-unit_volume.txt'
        out_path.write_text('1/m$^3$')

        out_path = cfg['OUT'] / 'multi_species-unit_surface.txt'
        out_path.write_text('1/m$^2$')

        _write_data_file(cfg, 'multi_species-values_volume.txt',
                         values=np.hstack((t_scale * t[:, np.newaxis], vol_seg_factor * N_vol)).T,
                         label=['t'] + [f'v{n + 1:02}' for n in range(c_count)], exp=True)

        _write_data_file(cfg, 'multi_species-values_surface.txt',
                         values=np.hstack((t_scale * t[:, np.newaxis], surf_seg_factor * N_surf)).T,
                         label=['t'] + [f's{n + 1:02}' for n in range(c_count)], exp=True)


def multi_species_bonds(cfg, t, N_surf_detail, N_vol_detail):
    """Export function for *incremental* and *fragmentation* multi species model simulations.

            :param cfg: simulation config
            :type cfg: dict
            :param t: fit time series
            :type t: ndarray
            :param N_surf_detail: number of molecules on the surface
            :type N_surf_detail: ndarray
            :param N_vol_detail: number of molecules in solution
            :type N_vol_detail: ndarray

            :result:
                * same files as :meth:`multi_species`
                * plot of the average excess bond count ``excess_bonds.pdf``
                * space separated data file containing the average number of excess bonds ``multi_species-excess_bonds``

            """
    c_count = N_surf_detail.shape[1]
    b_count = N_surf_detail.shape[2]
    line_cm = cm.get_cmap('plasma_r', c_count - 1)
    norm = matplotlib.colors.Normalize(vmin=1.5, vmax=c_count + 0.5)

    N_vol = np.sum(N_vol_detail, axis=2)
    N_surf = np.sum(N_surf_detail, axis=2)

    ab = ma.masked_array(np.zeros_like(N_vol))
    complete_detail = N_vol_detail + N_surf_detail
    all_excess_bonds = np.zeros(len(N_vol))
    for n_b in range(b_count):
        ab += n_b * complete_detail[:, :, n_b]
    average_excess_bonds = ab / (N_vol + N_surf)
    for n_c in range(c_count):
        all_excess_bonds += ab[:, n_c]
        average_excess_bonds[:, n_c][N_vol[:, n_c] < 0.1] = ma.masked

    # average excess bond
    f = plt.figure()
    ax1 = f.add_subplot(111)
    ax2 = ax1.twinx()
    t_scale, t_unit = data.Parameter.scale_time(np.max(t))
    ax1.set_ylabel(f"average number of excess bonds")
    ax2.set_ylabel(f"overall excess bonds (%)")
    ax1.set_xlabel(f"time $t$ ({t_unit})")
    for k in range(c_count-1):
        ax1.plot(t_scale * t, average_excess_bonds[:, k+1], c=line_cm(norm(k + 2)))

    ax2.plot(t_scale * t, all_excess_bonds/np.max(all_excess_bonds)*100, c='black')
    ax1.set_xlim(0, np.max(t_scale * t))
    ax1.set_ylim((0, ax1.get_ylim()[1]))
    ax2.set_ylim((0, ax1.get_ylim()[1]/(b_count-1)*100))
    cbar = f.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=line_cm), orientation='vertical',
                      label='carbon count', pad=0.11)
    cbar.ax.invert_yaxis()
    f.tight_layout()
    f.savefig(cfg['OUT'] / 'excess_bonds.pdf')
    plt.close(f)

    _write_data_file(cfg, 'multi_species-excess_bonds.txt',
                     values=np.hstack((t_scale * t[:, np.newaxis], ma.filled(average_excess_bonds, -1))).T,
                     label=['t'] + [f'b{n+1:02}' for n in range(c_count)])
    _write_data_file(cfg, 'multi_species-overall_excess_bonds.txt',
                     values=[t_scale * t, all_excess_bonds/np.max(all_excess_bonds)*100],
                     label=['t', 'oeb'])

    multi_species(cfg, t, N_surf, N_vol)