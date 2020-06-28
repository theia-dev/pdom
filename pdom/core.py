from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from scipy.optimize import minimize_scalar

from . import export
from .data import Parameter


class Simulate(object):
    """Class to simulation different degradation models

    :param config_file: .ini file to load
    :type config_file: str, Path
    :param out_folder: folder to save the results
    :type out_folder: str, Path, optional
    :param data_file: .json file containing experimental data
    :type data_file: str, Path, optional
    :param overwrites: overwrite settings from the config file
    :type overwrites: dict, optional
    :param resolution: time resolution for data export
    :type resolution: int, optional
    """
    def __init__(self, config_file, out_folder=None, data_file=None, overwrites=None, resolution=250):
        """Constructor method
        """
        self.base_dir = Path(__file__).absolute().parent
        self.parameter = Parameter(config_file=config_file, data_file=data_file, overwrites=overwrites)
        #: the simulation configuration created with :class:`pdom.data.Parameter` from the config file.
        self.cfg = self.parameter.config
        if out_folder is None:
            out_folder = Path('.').absolute() / self.cfg['SIMULATION']['id']
        out_folder_path = Path(out_folder)
        if out_folder_path.is_file():
            exit(f'ERROR: {out_folder_path.absolute()} is not a directory')
        out_folder_path.mkdir(parents=True, exist_ok=True)
        self.cfg['OUT'] = out_folder_path.absolute()
        self.resolution = resolution
        # simulation variables (filled in _prepare)
        #: time steps used by the solver - filled after :meth:`run`
        self.t = None
        #: time steps used for export - filled after :meth:`run`
        self.export_t = None

    def rhs_single(self, t, N, k):
        """Right-hand site for single species model.

        Calculate the derivative for a given concentration profile.

        :param N: number of molecules | 1-dim
        :type N: ndarray
        :param t: time (not used)
        :type t: float
        :param k: simulation constants (:math:`k_{\mathrm{ads}}`, :math:`k_{\mathrm{des}}`, :math:`k_{\mathrm{reac}}`)
        :type k: ndarray

        * ``N[0]`` is the number of molecules on the surface
        * ``N[1]`` is the number of molecules in solution

        :return: first derivative
        :rtype: ndarray

        """
        ka, kd, kr = k

        # Create empty array for the absolute flux per element
        F = np.zeros(2)

        # Calculate the surface concentrations
        C_AS = N[0] / self.cfg['CATALYST']['surface_total']

        # Calculate the volume concentration
        C_A = N[1] / self.cfg['CATALYST']['volume']
        theta = C_AS * self.cfg['MOLECULE']['molar_surface']

        # Calculate the fluxes
        J_ads = C_A * (1.0 - theta) * ka
        J_des = C_AS * kd
        J_reac = C_AS * kr

        F[0] = (-J_des + J_ads - J_reac) * self.cfg['CATALYST']['surface_total']
        F[1] = (+J_des - J_ads) * self.cfg['CATALYST']['surface_total']

        return F

    def rhs_multi(self, t, N_flat, k, N_shape):
        """Right-hand site for multi species models *fragmentation* and *incremental*.

        Calculate the derivative for a given concentration profile.

        :param N_flat: number of molecules | 2-dim flattened
        :type N_flat: ndarray
        :param t: time (not used)
        :type t: float
        :param k: simulation constants (:math:`k_{\mathrm{ads}}`, :math:`k_{\mathrm{des}}`, :math:`k_{\mathrm{reac}}`)
        :type k: ndarray
        :param N_shape: original shape of N
        :type N_shape: tuple

        * first dimension of N corresponds to the number of Carbon atoms ``carbon_count=index+1``
        * ``N[:, 0]`` is the number of molecules on the surface
        * ``N[:, 1]`` is the number of molecules in solution
        * the simulation constants must have the shape ``(3, max(carbon_count))``

        :return: first derivative
        :rtype: ndarray
        """
        N = N_flat.reshape(N_shape)
        ka, kd, kr = k

        # Create empty array for the absolute fluxes
        F = np.zeros(N_shape)

        # Calculate the concentrations
        C_S = N[:, 0] / self.cfg['CATALYST']['surface_total']
        C_L = N[:, 1] / self.cfg['CATALYST']['volume']
        theta = (sum(C_S[1:] * self.cfg['MOLECULE']['molar_surface_multi'][1:]) /
                 self.cfg['CATALYST']['surface_total'])

        if self.cfg['MULTI']['split_model'] == 'fragmentation':
            J_reac_in = np.zeros(N_shape[0])
            for i in range(1, N.shape[0]):
                for fragment_1 in range(int(np.floor((1 + i) / 2.0))):
                    fragment_2 = i - fragment_1 - 1
                    J_reac_in[fragment_1] += C_S[i] * kr[i] * 2.0 / i
                    if fragment_2 != fragment_1:
                        J_reac_in[fragment_2] += C_S[i] * kr[i] * 2.0 / i

        # Calculate the fluxes
        J_ads = C_L * (1 - theta) * ka
        J_des = C_S * kd
        J_reac_out = C_S * kr

        # Calculate the total fluxes
        F[:, 0] = (-J_des + J_ads) * self.cfg['CATALYST']['surface_total']
        F[:, 1] = (+J_des - J_ads) * self.cfg['CATALYST']['surface_total']

        F[0, :] = 0.0
        if self.cfg['MULTI']['split_model'] == 'incremental':
            F[:-1, 0] = F[:-1, 0] + J_reac_out[1:] * self.cfg['CATALYST']['surface_total']
            F[0, 0] = F[0, 0] + sum(J_reac_out[1:]) * self.cfg['CATALYST']['surface_total']
        elif self.cfg['MULTI']['split_model'] == 'fragmentation':
            F[:, 0] = F[:, 0] + J_reac_in[:] * self.cfg['CATALYST']['surface_total']

        F[1:, 0] = F[1:, 0] - J_reac_out[1:] * self.cfg['CATALYST']['surface_total']

        return F.flatten()

    def rhs_bonds(self, t, N_flat, k, N_shape):
        """Right-hand site for multi species model *excess_bonds*.

        Calculate the derivative for a given concentration profile.

        :param N_flat: number of molecules | 3-dim flattened
        :type N_flat: ndarray
        :param t: time (not used)
        :type t: float
        :param k: simulation constants (:math:`k_{\mathrm{ads}}`, :math:`k_{\mathrm{des}}`, :math:`k_{\mathrm{reac}}`)
        :type k: ndarray
        :param N_shape: original shape of N
        :type N_shape: tuple

        * first dimension of N corresponds to the number of Carbon atoms ``carbon_count=index+1``
        * second dimension of N corresponds to the number of excess bonds ``excess_bonds=index``
        * ``max(excess_bonds) < max(carbon_count)``
        * ``N[:, :, 0]`` is the number of molecules on the surface
        * ``N[:, :, 1]`` is the number of molecules in solution
        * the simulation constants must have the shape ``(3, max(excess_bonds)+1, max(carbon_count))``

        :return: first derivative
        :rtype: ndarray
        """

        N = N_flat.reshape(N_shape)

        ka, kd, kr = k

        # Create empty array for the absolute fluxes
        F = np.zeros(N_shape)

        # Calculate the concentrations
        C_S = N[:, :, 0] / self.cfg['CATALYST']['surface_total']
        C_L = N[:, :, 1] / self.cfg['CATALYST']['volume']
        theta = (sum(C_S[1:] * self.cfg['MOLECULE']['molar_surface_multi_eb'][1:]) /
                 self.cfg['CATALYST']['surface_total'])

        carbon_count = np.tile(range(N_shape[0]), (N_shape[1], 1)).T + 1
        excess_bonds = np.tile(range(N_shape[1]), (N_shape[0], 1))

        # Calculate the fluxes
        J_ads = C_L * (1 - theta) * ka
        J_des = C_S * kd
        J_reac_out = kr * C_S
        J_reac_in = np.zeros(N_shape[:-1])

        J_reac_in[1:, :-1] += (excess_bonds[1:, 1:])/(carbon_count[1:, 1:] - 1) * kr[1:, 1:] * C_S[1:, 1:]
        for n_c in range(2, N_shape[0]+1):
            for n_b in range(min(n_c, N_shape[1])):
                J_reac_in_part = 2 * (n_c - n_b - 1) / ((n_c - 1) ** 2) * kr[n_c - 1, n_b] * C_S[n_c - 1, n_b]
                if J_reac_in_part > 0:
                    for n_c_f1 in range(1, int(np.floor(n_c / 2))+1):
                        n_c_f2 = n_c - n_c_f1
                        n_b_f2 = max(int(np.round(n_c_f2 * n_b / n_c)), n_b-(n_c-n_c_f2-1))
                        n_b_f1 = n_b - n_b_f2

                        J_reac_in[n_c_f1 - 1, n_b_f1] += J_reac_in_part
                        if n_c_f1 != n_c_f2:
                            J_reac_in[n_c_f2 - 1, n_b_f2] += J_reac_in_part
                        pass

        # Calculate the total fluxes
        F[:, :, 0] = (-J_des + J_ads) * self.cfg['CATALYST']['surface_total']
        F[:, :, 1] = (+J_des - J_ads) * self.cfg['CATALYST']['surface_total']

        F[0, :, :] = 0.0

        F[:, :, 0] = F[:, :, 0] + J_reac_in[:] * self.cfg['CATALYST']['surface_total']
        F[1:, :, 0] = F[1:, :, 0] - J_reac_out[1:] * self.cfg['CATALYST']['surface_total']

        return F.flatten()

    def _calculate_error_interp(self, org_t, org_y, calc_t, calc_y):
        fit_y = np.interp(org_t, calc_t, calc_y)
        return self._calculate_error(org_y, fit_y)

    def _calculate_error(self, org_y, fit_y):
        if self.cfg['FIT']['search'] == "absolute":
            return np.sum(np.abs((fit_y - org_y)))
        elif self.cfg['FIT']['search'] == "relative":
            return np.sum(np.abs((fit_y - org_y) / fit_y))
        elif self.cfg['FIT']['search'] == "relative_square":
            return np.sum(((fit_y - org_y) / fit_y) ** 2)
        else:
            exit(f"FIT search type '{self.cfg['FIT']['search']}' does not exist.")

    def _prepare(self):
        self.t = [0, self.cfg['SIMULATION']['duration']]
        if self.resolution:
            self.export_t = np.linspace(0, self.cfg['SIMULATION']['duration'], self.resolution)
        else:
            self.export_t = np.linspace(0, self.cfg['SIMULATION']['duration'], 250)

    def _prepare_run(self, run_type):
        N_0_vol = self.cfg['SYSTEM']['concentration_solution']
        N_0_surf = None
        if 'concentration_surface' in self.cfg['SYSTEM']:
            N_0_surf = self.cfg['SYSTEM']['concentration_surface']
        if N_0_surf is None:
            N_0_surf = self.cfg['CATALYST']['surface_total'] / (
                    (self.cfg['SYSTEM']['k_des'] / self.cfg['SYSTEM']['k_ads']) * (
                        self.cfg['CATALYST']['volume'] / N_0_vol) +
                    self.cfg['MOLECULE']['molar_surface']
            )

        if run_type == 'single':
            k = [self.cfg['SYSTEM']['k_ads'],
                 self.cfg['SYSTEM']['k_des'],
                 self.cfg['SYSTEM']['k_reac']]
        else:
            k = [self.cfg['SYSTEM']['k_ads_multi'],
                 self.cfg['SYSTEM']['k_des_multi'],
                 self.cfg['SYSTEM']['k_reac_multi']]

        if run_type == 'bond':
            c_shape = (self.cfg['MOLECULE']['composition']["C"], self.cfg['MOLECULE']['excess_bonds'] + 1, 2)
            N_0 = np.zeros(c_shape)
            N_0[-1, -1, 0] = N_0_surf
            N_0[-1, -1, 1] = N_0_vol
        elif run_type == 'single':
            N_0 = np.array([N_0_surf, N_0_vol])
            c_shape = (2,)
        elif run_type == 'multi':
            c_shape = (self.cfg['MOLECULE']['composition']["C"], 2)
            N_0 = np.zeros(c_shape)
            N_0[-1, 0] = N_0_surf
            N_0[-1, 1] = N_0_vol
        else:
            raise KeyError(f'Unknown key {run_type} for run_type.')

        return N_0, k, c_shape

    def _run_single(self):
        print("Start calculating single species model")
        N_0, k, c_shape = self._prepare_run('single')

        result = solve_ivp(self.rhs_single, t_span=self.t, y0=N_0.flatten(), t_eval=self.export_t, args=(k,),
                           method='LSODA', atol=self.cfg['SOLVER']['atol'], rtol=self.cfg['SOLVER']['rtol'])
        c_vol, c_surf = result.y[1], result.y[0]

        export.single_species(self.cfg, self.export_t, c_surf, c_vol)
        print("Calculation finished!")
        print(f"Results saved in {self.cfg['OUT']}")

    def _run_multi(self):
        print("Start calculating multi species model")
        N_0, k, c_shape = self._prepare_run('multi')

        result = solve_ivp(self.rhs_multi, t_span=self.t, y0=N_0.flatten(),  t_eval=self.export_t, args=(k, N_0.shape),
                           method='LSODA', atol=self.cfg['SOLVER']['atol'], rtol=self.cfg['SOLVER']['rtol'])
        y_vec = result.y.T
        y_vec = y_vec.reshape(y_vec.shape[0], c_shape[0], c_shape[1])
        c_vol, c_surf = y_vec[:, :, 1], y_vec[:, :, 0]
        export.multi_species(self.cfg, self.export_t, c_surf, c_vol)
        print("Calculation finished!")
        print(f"Results saved in {self.cfg['OUT']}")

    def _run_bond(self):
        print("Start calculating multi species bond model")
        N_0, k, c_shape = self._prepare_run('bond')

        result = solve_ivp(self.rhs_bonds, t_span=self.t, y0=N_0.flatten(), t_eval=self.export_t, args=(k, N_0.shape),
                           method='LSODA', atol=self.cfg['SOLVER']['atol'], rtol=self.cfg['SOLVER']['rtol'])
        y_vec = result.y.T
        y_vec = y_vec.reshape(y_vec.shape[0], c_shape[0], c_shape[1], c_shape[2])
        c_vol, c_surf = y_vec[:, :, :, 1], y_vec[:, :, :, 0]
        export.multi_species_bonds(self.cfg, self.export_t, c_surf, c_vol)
        print("Calculation finished!")
        print(f"Results saved in {self.cfg['OUT']}")

    def _fit_dark_min(self, k_ads, final=False):
        sum_error = 0

        if final:
            t = self.export_t
            collect_results = []

        if self.cfg['DATA']['initial_concentration'].ndim == 2:
            N_vol_dark = self.cfg['DATA']['time_series'][:, 1, -1]
            N_surf_dark = np.sum(self.cfg['DATA']['initial_concentration'], axis=0)-self.cfg['DATA']['time_series'][:, 1, -1]
            k_des = k_ads * np.average(
                    N_vol_dark / self.cfg['CATALYST']['volume'] * (
                    self.cfg['CATALYST']['surface_total'] / N_surf_dark - self.cfg['MOLECULE']['molar_surface']
                )
            )
            k = (k_ads, k_des, 0)
            for n in range(self.cfg['DATA']['initial_concentration'].shape[1]):
                N_0 = np.array((self.cfg['DATA']['initial_concentration'][0, n],
                                self.cfg['DATA']['initial_concentration'][1, n]))
                if not final:
                    t = self.cfg['DATA']['time_series'][n][0]
                result = solve_ivp(self.rhs_single, t_span=self.t, y0=N_0.flatten(),
                                   t_eval=t, args=(k,),
                                   method='LSODA', atol=self.cfg['SOLVER']['atol'], rtol=self.cfg['SOLVER']['rtol'])
                y = result.y[1]
                if final:
                    collect_results.append(y)
                else:
                    sum_error += self._calculate_error(self.cfg['DATA']['time_series'][n][1], y)
        else:
            N_vol_dark = self.cfg['DATA']['time_series'][1, -1]
            N_surf_dark = self.cfg['DATA']['initial_concentration'][1] - self.cfg['DATA']['time_series'][1, -1]
            N_0 = np.array(self.cfg['DATA']['initial_concentration'])
            k_des = k_ads * (
                N_vol_dark/self.cfg['CATALYST']['volume']*(
                self.cfg['CATALYST']['surface_total']/N_surf_dark-self.cfg['MOLECULE']['molar_surface']
                )
            )
            k = (k_ads, k_des, 0)
            if not final:
                t = self.cfg['DATA']['time_series'][0]
            result = solve_ivp(self.rhs_single, t_span=self.t, y0=N_0.flatten(), t_eval=t, args=(k,),
                               method='LSODA', atol=self.cfg['SOLVER']['atol'], rtol=self.cfg['SOLVER']['rtol'])
            y = result.y[1]
            if final:
                collect_results.append(y)
            else:
                sum_error = self._calculate_error(self.cfg['DATA']['time_series'][1], y)

        if final:
            return self.export_t, collect_results, k_ads, k_des
        else:
            return sum_error

    def _fit_dark(self):
        print(f"Start fitting to data from dark experiment using '{self.cfg['FIT']['search']}' error estimation.")
        result_raw = minimize_scalar(self._fit_dark_min, bracket=(0, 1E-15), tol=1E-5, options={"maxiter": self.cfg['FIT']['max_step']})
        if result_raw.success:
            result = self._fit_dark_min(result_raw.x, final=True)
            export.fit_dark(self.cfg, result, result_raw)
            print(f"Fit finished after {result_raw.nit} iterations.")
            print(f"\tk_ads {result[2]:.3E} ")
            print(f"\tk_des {result[3]:.3E} ")
            print(f"\terror {result_raw.fun:.3E} ")
            print(f"Results saved in {self.cfg['OUT']}")
        else:
            print('ERROR: Fit was not successful.')
            print(result_raw)

    def _fit_single_min(self, k_reac, final=False):
        sum_error = 0

        k_ads = self.cfg['SYSTEM'].get('k_ads')
        k_des = self.cfg['SYSTEM'].get('k_des')

        if self.cfg['DATA']['time_series'].ndim == 3:
            N_vol_dark = np.average(self.cfg['DATA']['time_series'][:, 1, 0])
        else:
            N_vol_dark = self.cfg['DATA']['time_series'][1, 0]

        N_surf_dark = self.cfg['SYSTEM']['concentration_solution'] - N_vol_dark

        if k_des is None:
            k_des = k_ads * np.average(
                N_vol_dark / self.cfg['CATALYST']['volume'] * (
                        self.cfg['CATALYST']['surface_total'] / N_surf_dark - self.cfg['MOLECULE']['molar_surface']))

        if k_ads is None:
            k_ads = k_des / np.average(
                N_vol_dark / self.cfg['CATALYST']['volume'] * (
                        self.cfg['CATALYST']['surface_total'] / N_surf_dark - self.cfg['MOLECULE']['molar_surface']
                )
            )

        k = (k_ads, k_des, k_reac)

        N_0 = np.array((np.average(N_surf_dark),
                        np.average(N_vol_dark)))

        result = solve_ivp(self.rhs_single, t_span=self.t, y0=N_0.flatten(), t_eval=self.export_t, args=(k,),
                           method='LSODA', atol=self.cfg['SOLVER']['atol'], rtol=self.cfg['SOLVER']['rtol'])

        y = result.y[1]
        if final:
            collect_results = [y]
        else:
            if self.cfg['DATA']['time_series'].ndim == 3:
                for n in range(self.cfg['DATA']['time_series'].shape[0]):
                    sum_error += self._calculate_error_interp(self.cfg['DATA']['time_series'][n][0],
                                                              self.cfg['DATA']['time_series'][n][1],
                                                              self.export_t, y
                                                              )
            else:
                sum_error = self._calculate_error_interp(self.cfg['DATA']['time_series'][0],
                                                         self.cfg['DATA']['time_series'][1],
                                                         self.export_t, y
                                                         )

        if final:
            return self.export_t, collect_results, k_ads, k_des, k_reac
        else:
            return sum_error

    def _fit_single(self):
        print("Start fitting to data from reaction experiment.")
        result_raw = minimize_scalar(self._fit_single_min, bracket=(1E-8, 1E-4), tol=1E-5,
                                     options={"maxiter": self.cfg['FIT']['max_step']})
        if result_raw.success:
            result = self._fit_single_min(result_raw.x, final=True)
            export.fit_single(self.cfg, result, result_raw)
            print(f"Fit finished after {result_raw.nit} iterations.")
            print(f"\tk_ads: {result[2]:.3E} m/s")
            print(f"\tk_des: {result[3]:.3E} 1/s")
            print(f"\tk_reac: {result[4]:.3E} 1/s")
            print(f"\terror: {result_raw.fun:.3E} ")
            print(f"Results saved in {self.cfg['OUT']}")
        else:
            print('ERROR: Fit was not successful.')
            print(result_raw)

    def _fit_toc_func_reac(self, t, k_reac):
        return self._fit_toc_func(t, k_reac=k_reac)

    def _fit_toc_func_beta(self, t, beta_1):
        return self._fit_toc_func(t, beta_1=beta_1)

    def _fit_toc_func(self, t, k_reac=None, beta_1=None, full_resolution=False):
        if k_reac is not None:
            self.cfg['SYSTEM']['k_reac'] = k_reac
        if beta_1 is not None:
            self.cfg['MULTI_WEAK']['beta_1'] = beta_1
        try:
            del self.cfg['MULTI_WEAK']['beta_0']
        except KeyError:
            pass
        self.cfg = self.parameter.update_multi_k(self.cfg)

        if full_resolution:
            t = self.export_t

        if self.cfg['MULTI']['split_model'] == 'excess_bonds':
            N_0, k, c_shape = self._prepare_run('bond')

            result = solve_ivp(self.rhs_bonds, t_span=self.t, y0=N_0.flatten(), t_eval=t, args=(k, N_0.shape),
                               method='LSODA', atol=self.cfg['SOLVER']['atol'], rtol=self.cfg['SOLVER']['rtol'])
            y_vec = result.y.T
            y_vec = y_vec.reshape(y_vec.shape[0], c_shape[0], c_shape[1], c_shape[2])

            toc_sim = np.sum(y_vec[:, :, :, 1], axis=2)
            toc_sim = np.sum(toc_sim * (np.arange(c_shape[0]) + 1), axis=1)

        else:
            N_0, k, c_shape = self._prepare_run('multi')

            result = solve_ivp(self.rhs_multi, t_span=self.t, y0=N_0.flatten(), t_eval=t, args=(k, N_0.shape),
                               method='LSODA', atol=self.cfg['SOLVER']['atol'], rtol=self.cfg['SOLVER']['rtol'])
            y_vec = result.y.T
            y_vec = y_vec.reshape(y_vec.shape[0], c_shape[0], c_shape[1])

            toc_sim = np.sum(y_vec[:, :, 1] * (np.arange(c_shape[0]) + 1), axis=1)

        toc_sim = toc_sim/toc_sim[0]
        return toc_sim

    def _fit_toc(self):
        print("Start fitting to toc")
        t = self.cfg['DATA']['time_series'][0]
        toc = self.cfg['DATA']['time_series'][1] / self.cfg['DATA']['time_series'][1][0]

        if 'beta_1' in self.cfg['MULTI_WEAK']:
            popt, pcov = curve_fit(self._fit_toc_func_reac, t, toc, verbose=2, p0=0.01, bounds=(1E-6, 1), xtol=1.0E-4)
            fit_toc = self._fit_toc_func(t, popt[0], self.cfg['MULTI_WEAK']['beta_1'], full_resolution=True)

        elif 'k_reac' in self.cfg['SYSTEM']:
            popt, pcov = curve_fit(self._fit_toc_func_beta, t, toc, verbose=2, p0=0.1,  bounds=(1E-5, 10), xtol=1.0E-4)
            fit_toc = self._fit_toc_func(t, self.cfg['SYSTEM']['k_reac'], popt[0], full_resolution=True)

        else:
            raise KeyError('Either MULTI_WEAK|beta_1 or SYSTEM|k_reac needs to be set in the cunfiguration.')

        error_sd = np.sqrt(pcov[0, 0])

        export.fit_toc(self.cfg, error_sd, self.export_t, fit_toc)
        print(f"Fit finished")
        print(f"\tk_ads: {self.cfg['SYSTEM']['k_ads']:.3E} m/s")
        print(f"\tk_des: {self.cfg['SYSTEM']['k_des']:.3E} 1/s")
        print(f"\tk_reac: {self.cfg['SYSTEM']['k_reac']:.3E} 1/s")
        print(f"\tbeta_0: {self.cfg['MULTI_WEAK']['beta_0']:.3E} 1/s")
        print(f"\tbeta_1: {self.cfg['MULTI_WEAK']['beta_1']:.3E} 1/s")
        print(f"\terror: {error_sd:.3E} ")
        print(f"Results saved in {self.cfg['OUT']}")

    def run(self):
        """ Runs either a simulation or performs a parameter
        fit based on the information stored in :attr:`cfg`.
        The results are saved to the file system.
        """
        self._prepare()
        if self.cfg['SIMULATION']['fit']:
            if self.cfg['FIT']['type'] == 'dark':
                self._fit_dark()
            elif self.cfg['FIT']['type'] == 'reac':
                self._fit_single()
            elif 'toc' in self.cfg['FIT']['type']:
                self._fit_toc()
        else:
            if self.cfg['SIMULATION']['multi']:
                if self.cfg['MULTI']['split_model'] == 'excess_bonds':
                    self._run_bond()
                else:
                    self._run_multi()
            else:
                self._run_single()
