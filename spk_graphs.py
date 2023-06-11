import os
import numpy as np
import random
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.legend import Legend
from matplotlib.legend_handler import HandlerBase, HandlerPatch, update_from_first_child
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import json


def spk_graphs(df_out, file_name="", dict_colors_groups=None,
               dict_colors_intergroup_stars=None, autocomplete_colors=False,
               independent=False, gf_uw2="circumferences", asso_par=None):
    """
    Function to plot the Boettlinger diagrams (the UV and UW planes), with and
    without zoom in on the Eggen's young disk, and the Toomre diagram.

    Parameters
    ----------

    df_out :  DataFrame
        The DataFrame resulting from spk_groups or, at least, one containing
        the following columns: U, V, and W.
    file_name : str, optional
        A string to be included in the names of the output figures. The default
        value is empty.
    dict_colors_groups : dict or dict-like str, optional
    dict_colors_intergroup_stars : dict or dict-like str, optional
        dict_colors_groups and dict_colors_intergroup_stars are python
        dictionaries or dictionary-like strings containing the name (key) and
        color (value) given by the user to the stars that fall in one
        association/group or between several associations/groups apart from
        those provided by default. To "remove" a default
        association/group/intergroup from the graphics, the corresponding
        dictionary must contain "default name": None (for example, to avoid
        showing LA stars in red, dict_colors_groups must contain "LA": None).
    autocomplete_colors : boolean, optional
        In case the input DataFrame contains stars belonging to any new
        association, group, or intergroup and it has not been defined in the
        corresponding dictionary, the "autocomplete_colors" parameter sets
        whether showing these stars in a randomly generated color (True) or
        displaying them according to their stellar population if given (False).
        Caution: autocomplete_colors only avoid repetition if the colors are
        given by their hexadecimal name.
    independent : boolean, optional
        Different figures for each UV and WV plane? Default is False.
    gf_uw2 : str, optional
        The gf_uw2 parameter sets the aspect ratio of the Toomre diagram. The
        values, circumferences (default) or ellipses, indicate how lines of
        constant total velocity look in the resulting diagram.
    asso_par : str, optional
        Name (without extension) of the custom file for the parameters of the
        associations and moving groups. The CSV file must have the same
        structure as the default provided and be located in the current
        directory.
    """

    cols = {'Vr_kms-1', 'eVr_kms-1', 'U_kms-1', 'eU_kms-1', 'V_kms-1',
            'eV_kms-1', 'W_kms-1', 'eW_kms-1'}
    for col in cols:
        df_out[col] = df_out[col].astype(float)

    # CSV file containing the parameters of groups and associations
    if asso_par is None:
        try:
            df_par = pd.read_csv(
                os.path.join(os.path.dirname(__file__),
                             'association_parameters.csv'), index_col=None)
        except FileNotFoundError as e:
            raise FileNotFoundError('The file association_parameters.csv ' +
                                    'was not found in the corresponding ' +
                                    'directory.').with_traceback(e.__traceback__)
    else:
        try:
            df_par = pd.read_csv(asso_par+'.csv', index_col=None)
        except FileNotFoundError as e:
            raise FileNotFoundError(
                'The file '+ asso_par +'.csv was not found in the current'+
                ' directory.').with_traceback(e.__traceback__)

    # Parse str to dict if necessary
    class InvalidDictLikeStr(Exception):
        """Exception raised for errors in the input dict-like strings.

        Attributes:
            message -- explanation of the error
        """

        def __init__(self, message):
            super().__init__(message)

    def parseDict(strDict, error, clss):
        try:
            json.loads(strDict)
        except ValueError:
            raise InvalidDictLikeStr(error)
            # return None
        if clss == 'group':
            if any(('/' in key) for key in json.loads(strDict)) is True:
                raise InvalidDictLikeStr(error)
                # return None
        if clss == 'intergroup':
            if any(('/' in key) for key in json.loads(strDict)) is False:
                raise InvalidDictLikeStr(error)
                # return None
        d_out = json.loads(strDict)
        for key, value in d_out.items():
            if value == "None":
                d_out[key] = None
        return d_out

    if type(dict_colors_groups) == str:
        errorg = str(" Invalid dictionary-like string: dict_colors_groups" +
                     " is a python dictionary or a dictionary-like string in" +
                     " the form '{\"group\": \"Value\"}', where value is" +
                     " either a hexadecimal color or None (for removing" +
                     " default colors)")
        dict_colors_groups = parseDict(dict_colors_groups, errorg, 'group')

    if type(dict_colors_intergroup_stars) == str:
        errori = str("Invalid dictionary-like string:" +
                     " dict_colors_intergroup_stars is a python dictionary" +
                     " or a dictionary-like string in the form" +
                     " '{\"intergroup\": \"Value\"}', where intergroup is in" +
                     " the form \"group1/group2(/group3..)\" and value is" +
                     " either a hexadecimal color or None (for removing" +
                     " default colors)")
        dict_colors_intergroup_stars = parseDict(dict_colors_intergroup_stars,
                                                 errori, 'intergroup')

    # Default dicts for colors
    dict_colors_groups_default = {'LA': '#ff0000', 'CAS': '#00ffff',
                                  'IC': '#ffc0cb', 'HS': '#000080',
                                  'UMA': '#006400'}

    dict_colors_intergroup_stars_default = {'LA/IC': '#db4bda',
                                            'LA/CAS': '#aaa662',
                                            'LA/IC/CAS': '#c9ff27',
                                            'IC/CAS': '#3e82fc'}

    # Removing default associations, groups or intergroups and appending new
    # ones.
    def edit_colors(d1, d2):
        if d2 != None:
            dout = {**d1, **d2}
            keys_to_remove = []
            for k, v in dout.items():
                if v is None:
                    keys_to_remove.append(k)
            for key in keys_to_remove:
                del dout[key]
        else:
            dout = d1
        return dout

    dict_col_groups = edit_colors(dict_colors_groups_default,
                                  dict_colors_groups)
    dict_col_inter = edit_colors(dict_colors_intergroup_stars_default,
                                 dict_colors_intergroup_stars)

    if autocomplete_colors is True:
        col_list = []
        new_groups = []
        col_list_inter = []
        new_intergroups = []

        def random_color():
            col = "#"+''.join([random.choice(
                              '0123456789ABCDEF') for j in range(6)])
            while (col in dict_col_groups.items() or
                   col in dict_col_inter.items() or
                   col in col_list or
                   col in col_list_inter or
                   col == '#000000' or col == '#808080' or
                   col == '#FFA500'):
                col = "#"+''.join([random.choice(
                    '0123456789ABCDEF') for j in range(6)])
            return col

        for key in df_out.groupby('Kin_Pop').groups:
            if (key not in dict_col_groups.keys() and
                key not in dict_col_inter.keys() and
                key not in dict_colors_groups_default.keys() and
                key not in dict_colors_intergroup_stars_default.keys() and
                    key != 'NNN'):
                col = random_color()
                if '/' not in key:
                    new_groups.append(key)
                    col_list.append(str(col))
                else:
                    new_intergroups.append(key)
                    col_list_inter.append(str(col))

        if new_groups != []:
            dict_col_new_groups = {}
            for k, v in zip(new_groups, col_list):
                dict_col_new_groups[k] = v
            # In case a new group has been defined, but any star falls into
            # it.
            for i in range(len(df_par['ASSOCIATION/GROUP'])):
                if (df_par['ASSOCIATION/GROUP'][i] not in
                    dict_col_groups.keys() and
                    df_par['ASSOCIATION/GROUP'][i] not in
                        dict_colors_groups_default.keys()):
                    col = random_color()
                    dict_col_new_groups[
                        df_par['ASSOCIATION/GROUP'][i]] = col
            dict_col_groups = {**dict_col_groups, **dict_col_new_groups}

        if new_intergroups != []:
            dict_col_new_inter = {}
            for k, v in zip(new_intergroups, col_list_inter):
                dict_col_new_inter[k] = v
            dict_col_inter = {**dict_col_inter, **dict_col_new_inter}

    # Function to plot the limits of the Eggen's definition of the young disk
    def eggens_yd(subplot_ax, plane):
        ax = subplot_ax
        if plane == 'UV':
            u_1 = np.arange(-48, -7, 0.05)
            ax.plot(u_1, -30*(u_1/u_1), linestyle='--', color='#000000',
                    linewidth=0.7, zorder=14)
            u_2 = np.arange(-7, 20, 0.05)
            ax.plot(u_2, (u_2 + 7.0*(u_2/u_2))*((-5.0+30.0)/(20.0+7.0)) - 30.0,
                    linestyle='--', color='#000000', linewidth=0.7, zorder=14)
            v_3 = np.arange(-5.0, 2.7, 0.05)
            ax.plot(20.0*(v_3/v_3), v_3, linestyle='--', color='#000000',
                    linewidth=0.7, zorder=14)
            u_4 = np.arange(8., 20., 0.05)
            ax.plot(u_4, 2.7*(u_4/u_4), linestyle='--', color='#000000',
                    linewidth=0.7, zorder=14)
            u_5 = np.arange(0., 8., 0.05)
            ax.plot(u_5, (1/8) * u_5 * (2.7+2.0) - 2, linestyle='--',
                    color='#000000', linewidth=0.7, zorder=14)
            u_6 = np.arange(0., -28., -0.05)
            ax.plot(u_6, -2.0*np.ones(len(u_6)), linestyle='--',
                    color='#000000', linewidth=0.7, zorder=14)
            u_7 = np.arange(-28., -48., -0.05)
            ax.plot(u_7, -2 - (u_7+28.0) * ((-2.0+14.5)/(-48.5+28.0)),
                    linestyle='--', color='#000000', linewidth=0.7, zorder=14)
            v_8 = np.arange(-2 - (-48.0+28.0) * ((-2.0+14.5)/(-48.5+28.0)),
                            -30., -0.05)
            ax.plot(-48.0 * (v_8/v_8), v_8, linestyle='--', color='#000000',
                    linewidth=0.7, zorder=14)
        elif plane == 'WV':
            w_1 = np.arange(-30, 15, 0.05)
            ax.plot(w_1, -30*(w_1/w_1), linestyle='--', color='k',
                    linewidth=0.7, zorder=14)
            ax.plot(w_1, 2*(w_1/w_1),
                    linestyle='--', color='k', linewidth=0.7, zorder=14)
            v_4 = np.arange(-30.0, 2.0, 0.05)
            ax.plot(-30*(v_4/v_4), v_4, linestyle='--', color='k',
                    linewidth=0.7, zorder=14)
            ax.plot(15*(v_4/v_4), v_4, linestyle='--', color='k',
                    linewidth=0.7, zorder=14)
        return ax

    # Fuction to plot each plane of the Boettlinger diagrams
    def new_ax(ax, x, y, df_out, df_par, dict_col_groups,
               dict_col_inter, autocomplete_colors, sizes, xlabel, x1, x2,
               ylabel, y1, y2):
        # Stellar populations
        if ('Kin_Pop' in df_out.columns.values) is True:
            # Thin disk stars (D)- open grey circles
            ax.scatter(df_out[df_out['Kin_Pop'] == "D"][x],
                       df_out[df_out['Kin_Pop'] == "D"][y], s=sizes['D'],
                       facecolors='none', edgecolors='#808080', linewidth=0.25,
                       zorder=1)
            # Stars between the thin and thick disks (D-TD): grey squares with
            # black borders
            ax.scatter(df_out[df_out['Kin_Pop'] == "TD-D"][x],
                       df_out[df_out['Kin_Pop'] == "TD-D"][y], s=sizes['D-TD'],
                       marker='s', facecolors='#808080', edgecolors='black',
                       linewidth=0.25, zorder=14)
            # Thick disk stars (TD)- black filled circles
            ax.scatter(df_out[df_out['Kin_Pop'] == "TD"][x],
                       df_out[df_out['Kin_Pop'] == "TD"][y], s=sizes['TD'],
                       facecolors='#000000', edgecolors='#000000', zorder=14)
            # Halo stars (H)- black stars
            ax.scatter(df_out[df_out['Kin_Pop'] == "H"][x],
                       df_out[df_out['Kin_Pop'] == "H"][y], marker='*',
                       s=sizes['H'], color='#000000', linewidth=0.25,
                       zorder=14)

        # YD stars- open orange circles (caution: if the star belongs to a
        # group or association is displayed using the corresponding color of
        # that group or association)
        if ('YD' in df_out.columns.values) is True:
            ax.scatter(df_out[df_out['YD'] == "YD"][x],
                       df_out[df_out['YD'] == "YD"][y], s=sizes['YD'],
                       facecolors='none', edgecolors='#FFA500', linewidth=0.25,
                       zorder=1)

        # Candidate members of moving groups and associations and intergroup
        # stars
        if ('Kin_Group' in df_out.columns.values) is True:
            t = 0
            for key, value in dict_col_groups.items():
                ax.scatter(df_out[df_out['Kin_Group'] == key][x],
                           df_out[df_out['Kin_Group'] == key][y],
                           s=sizes['groups'], facecolors='none',
                           edgecolors=value, zorder=2+t)
                t += 1

            t = 0
            for key, value in dict_col_inter.items():
                ax.scatter(df_out[df_out['Kin_Group'] == key][x],
                           df_out[df_out['Kin_Group'] == key][y],
                           s=sizes['inter'], facecolors='none',
                           edgecolors=value, zorder=2+t)
                t += 1

        # Centers of the moving groups and associations: filled crosses
        for i in range(len(df_par)):
            if df_par['ASSOCIATION/GROUP'][i] in dict_col_groups.keys():
                ax.scatter(df_par[x][i], df_par[y][i], marker='+',
                           s=sizes['cross'], facecolors=dict_col_groups[
                    df_par['ASSOCIATION/GROUP'][i]], zorder=7+i)

        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(which='both', direction='in', top=True, bottom=True,
                       left=True, right=True)

        ax.set_xlim(x1, x2)
        ax.set_ylim(y1, y2)
        return ax

    params = {'lines.markersize': 1.5}
    plt.rcParams.update(params)
    scale = plt.rcParams['lines.markersize']**2

    # Dicts for the size of each marker in the new_ax function
    sizes_zoom = {'D': 7.*scale, 'D-TD': 12.*scale, 'TD': 8.*scale,
                  'H': 60.*scale, 'YD': 7.*scale, 'groups': 7.*scale,
                  'inter': 7.*scale, 'cross': 500.*scale}
    sizes_full = {'D': 7.*scale, 'D-TD': 12.*scale, 'TD': 8.*scale,
                  'H': 60.*scale, 'YD': 7.*scale, 'groups': 7.*scale,
                  'inter': 7.*scale, 'cross': 150.*scale}

    # Boettlinger diagrams
    if file_name is not "":
        out_name = 'SteParKin_' + file_name + '_'
    else:
        out_name = 'SteParKin_'

    if independent is False:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11., 11.))

    # V vs U zoomed in on the YD
    if independent is True:
        fig1, ax1 = plt.subplots(1, 1, figsize=(5.5, 5.5))

    ax1 = new_ax(ax1, 'U_kms-1', 'V_kms-1', df_out, df_par, dict_col_groups,
                 dict_col_inter, autocomplete_colors, sizes_zoom, 'U [km/s]',
                 -55.0, 30.0, 'V [km/s]', -35.0, 12.0)
    ax1 = eggens_yd(ax1, 'UV')

    if independent is True:
        name_fig1 = out_name + 'UV_zoom.pdf'
        fig1.savefig(name_fig1, bbox_inches='tight', dpi=600)

    # V vs W zoomed in on the YD
    if independent is True:
        fig2, ax2 = plt.subplots(1, 1, figsize=(5.5, 5.5))

    ax2 = new_ax(ax2, 'W_kms-1', 'V_kms-1', df_out, df_par, dict_col_groups,
                 dict_col_inter, autocomplete_colors, sizes_zoom, 'W [km/s]',
                 -50.0, 20.0, 'V [km/s]', -50.0, 15.0)
    ax2 = eggens_yd(ax2, 'WV')
    if independent is True:
        name_fig2 = out_name + 'WV_zoom.pdf'
        fig2.savefig(name_fig2, bbox_inches='tight', dpi=600)

    # V vs U (no zoom)
    u_min = np.floor(min(df_out['U_kms-1'])) - 10.
    if u_min > -58.0:
        u_min = -58.0
    u_max = np.ceil(max(df_out['U_kms-1'])) + 10.
    if u_max < 30.0:
        u_max = 30.0
    v_min = np.floor(min(df_out['V_kms-1'])) - 10.
    if v_min > -40.0:
        v_min = -40.0
    v_max = np.ceil(max(df_out['V_kms-1'])) + 10.
    if v_max < 15.0:
        v_max = 15.0

    if independent is True:
        fig3, ax3 = plt.subplots(1, 1, figsize=(5.5, 5.5))

    ax3 = new_ax(ax3, 'U_kms-1', 'V_kms-1', df_out, df_par, dict_col_groups,
                 dict_col_inter, autocomplete_colors, sizes_full, 'U [km/s]',
                 u_min, u_max, 'V [km/s]', v_min, v_max)
    ax3 = eggens_yd(ax3, 'UV')
    
    if independent is True:
        name_fig3 = out_name + 'UV.pdf'
        fig3.savefig(name_fig3, bbox_inches='tight', dpi=600)

    # V vs. W (no zoom)
    w_min = np.floor(min(df_out['W_kms-1'])) - 10.
    if w_min > -58.0:
        w_min = -58.0
    w_max = np.ceil(max(df_out['W_kms-1'])) + 10.
    if w_max < 30.0:
        w_max = 30.0
    v_min = np.floor(min(df_out['V_kms-1'])) - 10.
    if v_min > -40.0:
        v_min = -40.0
    v_max = np.ceil(max(df_out['V_kms-1'])) + 10.
    if v_max < 15.0:
        v_max = 15.0

    if independent is True:
        fig4, ax4 = plt.subplots(1, 1, figsize=(5.5, 5.5))

    ax4 = new_ax(ax4, 'W_kms-1', 'V_kms-1', df_out, df_par, dict_col_groups,
                 dict_col_inter, autocomplete_colors, sizes_full, 'W [km/s]',
                 w_min, w_max, 'V [km/s]', v_min, v_max)
    ax4 = eggens_yd(ax4, 'WV')
    if independent is True:
        name_fig4 = out_name + 'WV.pdf'
        fig4.savefig(name_fig4, bbox_inches='tight', dpi=600)

    if independent is False:
        name_fig1 = out_name + 'UV_WV.pdf'
        fig.savefig(name_fig1, bbox_inches='tight', dpi=600)

    # Correction to a dynamical local standard of rest (LSR):
    # The Suns peculiar motion relative to the LSR is assumed to be
    # (U, V, W) = (+10.00,+5.25,+7.17) km s-1 (Dehnen & Binney 1998)
    u_s = +10.00
    v_s = +5.25
    w_s = +7.17
    u_LSR = df_out['U_kms-1'] + u_s*(np.ones(len(df_out['U_kms-1'])))
    v_LSR = df_out['V_kms-1'] + v_s*(np.ones(len(df_out['V_kms-1'])))
    w_LSR = df_out['W_kms-1'] + w_s*(np.ones(len(df_out['W_kms-1'])))
    uw2 = (u_LSR**(2)+w_LSR**(2))**(0.5)

    dict2 = {'u_LSR': u_LSR, 'v_LSR': v_LSR, 'w_LSR': w_LSR, 'uw2': uw2}
    df2 = pd.DataFrame(dict2)

    # Toomre Diagram: V vs. uw2
    if max(v_LSR - 10*np.ones(len(v_LSR))) > 0.:
        x_max = max(v_LSR + 10*np.ones(len(v_LSR)))
    else:
        x_max = 10.

    x_min = min(v_LSR - 10*np.ones(len(v_LSR)))
    if (df2['uw2'].max() - df2['uw2'].min() < 50.):
        y_min = -0.1
    else:
        y_min = -5.
    y_max = max(uw2 + 10**np.ones(len(uw2)))

    # Aspect ratio of the figure depends on the option chosen.
    # Vt=cte lines look like circumferences (default) or ellipses.
    if gf_uw2 == 'ellipses':
        he = 5.5
        w = 5.5
    elif gf_uw2 == 'circumferences':
        if (x_max-x_min) > (y_max-y_min):
            aspect_ratio = (y_max-y_min)/(x_max-x_min)
            w = 5.5
            he = w * aspect_ratio
        else:
            aspect_ratio = (x_max-x_min)/(y_max-y_min)
            he = 5.5
            w = he * aspect_ratio
    else:
        raise ValueError('The options for gf_uw2 are ellipses or'
                         + ' circumferences.')

    fig2, ax1 = plt.subplots(1, 1, figsize=(w, he))
    # Stellar populations
    if ('Kin_Pop' in df_out.columns.values) is True:
        # D
        ax1.scatter(df2[df_out['Kin_Pop'] == "D"]['v_LSR'],
                    df2[df_out['Kin_Pop'] == "D"]['uw2'], s=5.*scale,
                    facecolors='none', edgecolors='#808080', linewidth=0.25,
                    zorder=1)
        # D-TD
        ax1.scatter(df2[df_out['Kin_Pop'] == "TD-D"]['v_LSR'],
                    df2[df_out['Kin_Pop'] == "TD-D"]['uw2'], s=5.*scale,
                    marker='s', facecolors='#808080', edgecolors='black',
                    linewidth=0.25, zorder=1)
        # TD
        ax1.scatter(df2[df_out['Kin_Pop'] == "TD"]['v_LSR'],
                    df2[df_out['Kin_Pop'] == "TD"]['uw2'], s=2.5*scale,
                    facecolors='#000000', edgecolors='#000000', zorder=1)
        # H
        ax1.scatter(df2[df_out['Kin_Pop'] == "H"]['v_LSR'],
                    df2[df_out['Kin_Pop'] == "H"]['uw2'], marker='*',
                    s=60.*scale, color='#000000', linewidth=0.25, zorder=1)

    # YD stars
    if ('YD' in df_out.columns.values) is True:
        ax1.scatter(df2[df_out['YD'] == "YD"]['v_LSR'],
                    df2[df_out['YD'] == "YD"]['uw2'], s=5.*scale,
                    facecolors='none', edgecolors='#FFA500', linewidth=0.25,
                    zorder=1)

    # Centers of moving groups and associations
    for i in range(len(df_par)):
        if df_par['ASSOCIATION/GROUP'][i] in dict_col_groups.keys():
            ax1.scatter(df_par['V_kms-1'][i] + v_s,
                        ((df_par['U_kms-1'][i]+u_s)**(2)
                         + (df_par['W_kms-1'][i]+w_s)**(2))**(0.5), marker='+',
                        s=125.0*scale, facecolors=dict_col_groups[
                df_par['ASSOCIATION/GROUP'][i]],
                zorder=7+i)

    # Candidate members of moving groups and associations and intergroup stars
    if ('Kin_Group' in df_out.columns.values) is True:
        t = 0
        for key, value in dict_col_groups.items():
            ax1.scatter(df2[df_out['Kin_Group'] == key]['v_LSR'],
                        df2[df_out['Kin_Group'] == key]['uw2'], s=5.*scale,
                        facecolors='none', edgecolors=value, linewidth=0.25,
                        zorder=2+t)
            t += 1

        t = 0
        for key, value in dict_col_inter.items():
            ax1.scatter(df2[df_out['Kin_Group'] == key]['v_LSR'],
                        df2[df_out['Kin_Group'] == key]['uw2'], s=5.*scale,
                        facecolors='none', edgecolors=value, linewidth=0.25,
                        zorder=2+t)
            t += 1

    # Lines of constant total velocity (Vt= cte)
    def vt_circles(aa, bb):
        X0 = 0.0
        Y0 = 0.0
        XR = np.zeros(360)
        YR = np.zeros(360)
        for i in range(len(XR)):
            RA = i*(np.pi/180.0)
            RR = (aa**(2))*(bb**(2))
            RR = RR/((aa**(2))*(np.sin(RA))**(2)
                     + (bb**(2))*(np.cos(RA))**(2))
            RR = RR**(0.5)
            XR[i] = X0 + RR*np.cos(RA)
            YR[i] = Y0 + RR*np.sin(RA)
        return XR, YR

    aa = 25.0
    bb = 25.0
    ax1.plot(vt_circles(aa, bb)[0], vt_circles(aa, bb)[1], linestyle='--',
             color='#000000', linewidth=0.5, zorder=6)
    aaa = np.arange(50.0, max(df2['u_LSR']) + 100.0, 50.0)
    bbb = np.arange(50.0, max(df2['uw2']) + 100.0, 50.0)
    for i in range(len(aaa)):
        ax1.plot(vt_circles(aaa[i], bbb[i])[0], vt_circles(aaa[i], bbb[i])[1],
                 linestyle='--', color='#000000', linewidth=0.5, zorder=6)

    ax1.set_ylabel('(U$^2$+W$^2$)$^{1/2}$ [km/s]', fontsize=9)
    ax1.set_xlabel('V [km/s]', fontsize=9)
    ax1.tick_params(which='both', direction='in')
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    name_fig2 = out_name + 'UW2_V.pdf'
    fig2.savefig(name_fig2, bbox_inches='tight', dpi=600)

    # SteParKin legend
    yd = Line2D([], [], color='#FFA500', marker='o', linestyle='None',
                linewidth=0.25, markersize=7., markerfacecolor='none',
                label=r"$\tt{Eggen's \ YD}$")
    d = Line2D([], [], color='#808080', marker='o', linestyle='None',
               linewidth=0.25, markersize=7., markerfacecolor='none',
               label='D')
    tdd = Line2D([], [], color='#000000', marker='s', linestyle='None',
                 linewidth=0.25, markersize=12., markerfacecolor='#808080',
                 label='TD-D')
    td = Line2D([], [], color='#000000', marker='.', linestyle='None',
                linewidth=0.25, markersize=8., markerfacecolor='#000000',
                label='TD')
    h = Line2D([], [], color='#000000', marker='*', linestyle='None',
               linewidth=0.25, markersize=12., markerfacecolor='#000000',
               label='H')

    if dict_col_groups is not None:
        groups_legend = []
        for key, value in dict_col_groups.items():
            groups_legend.append(Line2D([], [], color=value, marker='o',
                                        linestyle='None', linewidth=0.25,
                                        markersize=7., markerfacecolor='none',
                                        label=key))
    if dict_col_inter is not None:
        inter_legend = []
        for key, value in dict_col_inter.items():
            inter_legend.append(Line2D([], [], color=value, marker='o',
                                       linestyle='None', linewidth=0.25,
                                       markersize=7., markerfacecolor='none',
                                       label=key))

    fig_leg, ax = plt.subplots(figsize=(2, 0.75))

    legend_handles = []
    sub_ti = []

    if ('Kin_Pop' in df_out.columns.values) is True:
        title_pop = Line2D([], [], linewidth=0,
                           label=r"$\tt{Kin. Populations}$")
        legend_handles += [title_pop]
        populations_legend = [d, tdd, td, h]
        sub_ti += [title_pop]
        legend_handles += populations_legend
    if ('YD' in df_out.columns.values) is True:
        yd_legend = [yd]
        legend_handles += yd_legend
    if groups_legend is not None:
        title_g = Line2D([], [], linewidth=0,
                         label=r"$\tt{Kin. Groups}$")
        legend_handles += [title_g]
        sub_ti += [title_g]
        legend_handles += groups_legend
    if inter_legend is not None:
        title_i = Line2D([], [], linewidth=0,
                         label=r"$\tt{Intergroups}$")
        legend_handles += [title_i]
        sub_ti += [title_i]
        legend_handles += inter_legend

    ncol = 1
    if 10 < len(legend_handles) <= 15:
        ncol = 2
    elif 15 < len(legend_handles) <= 20:
        ncol = 3
    else:
        ncol = 4

    # The following treatment of the subtitles is from
    # https://github.com/mwaskom/seaborn/issues/1440
    def subtitle_handler_factory(inherit_from):
        """Class factory to subclass Handlers and add our custom functionality
        """
        class SubtitleHandler(inherit_from):
            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                handlebox.set_visible(False)
                return inherit_from.legend_artist(self, legend,
                                                  orig_handle, fontsize,
                                                  handlebox)

        # HandlerPatch class needs a special unpdate_func
        if inherit_from is HandlerPatch:
            return SubtitleHandler(update_func=update_from_first_child)
        return SubtitleHandler()

    def subtitle_handler_map(subtitles):
        defaults_handler_map = Legend.get_default_handler_map()
        handler_map = {}

        for orig_handle in subtitles:
            handler = Legend.get_legend_handler(
                defaults_handler_map, orig_handle)

            # Subclass the Handler
            new_handler = subtitle_handler_factory(type(handler))
            handler_map[orig_handle] = new_handler
        return handler_map

    handler_map = subtitle_handler_map(sub_ti)

    legend = ax.legend(handles=legend_handles, loc='lower left',
                       handler_map=handler_map, bbox_to_anchor=(-0.5, -0.75),
                       ncol=ncol, title=r"$\tt{\bf SteParKin}$", )
    legend.get_title().set_fontsize('12')
    ax.axis('off')
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
                        hspace=0, wspace=0)
    plt.margins(0, 0)
    fig_leg.savefig(out_name+'leg.png', pad_inches=0.1,
                    bbox_inches='tight', dpi=600)
