from spk_groups import spk_groups
from spk_graphs import spk_graphs
import pandas as pd
import numpy as np
import argparse
import ast


def spk_wrapper(input_file, errors="drop", file_name="",
                dict_colors_groups=None, dict_colors_intergroup_stars=None,
                autocomplete_colors=False, independent=False,
                gf_uw2="circumferences", asso_par=None):
    """
    This function:
        # Reads the input file (.csv or .txt)
        # Checks if the required columns are provided
        # Calls spk_groups (calculations) and spk_graphs (plotting graphs)

    Parameters
    ----------

    input_file : str
        Name (with extension) of the input file. The extension must be .txt or
        .csv. This file must contain the following columns:
            Name : str
                Name of the star
            RA_J2000_deg, DE_J2000_deg : float
                Equatorial coordinates (degrees) of the star in J2000 and
                epoch 2000
            muRA_masa-1, emuRA_masa-1, muDE_masa-1, emuDE_masa-1 : float
                Proper motions and their errors in mas yr-1
            Vr_kms-1, eVr_kms-1 : float
                Radial velocity and its error in km s-1
        and one of this two options:
            1) pi_mas and epi_mas : Parallax and its error in mas
            2) d_pc and ed_pc : Distance and its error in pc
        The names of all the columns are case sensitive. Not any column order
        is required.
    errors : str
        This parameter sets the behavior when an error is zero or NA/NAN.
        Default is "drop": the function skips any star with any error equals
        zero or NA/NAN (no record is left in the output file). Although,
        strictly speaking, only EPLX should not be zero to avoid division by
        zero when calculating the errors of the Galactic-space velocity
        components, we set this option as default because those errors may be
        artificially low otherwise. "zeros" keeps the stars with errors equals
        zero or NA/NAN  but converts NA/NAN values into zeros.
    file_name : str
        A string to be included in the names of the output figures. The default
        value is empty.
    dict_colors_groups : dict or dict-like str
    dict_colors_intergroup_stars : dict or dict-like str
        dict_colors_groups and dict_colors_intergroup_stars are python
        dictionaries or dictionary-like strings containing the name (key) and
        color (value) given by the user to the stars that fall in one
        association/group or between several associations/groups apart from
        those provided by default. To "remove" a default
        association/group/intergroup from the graphics, the corresponding
        dictionary must contain "default name": None (for example, to avoid
        showing LA stars in red, dict_colors_groups must contain "LA": None).
    autocomplete_colors : boolean
        In case the input DataFrame contains stars belonging to any new
        association, group, or intergroup and it has not been defined in the
        corresponding dictionary, the "autocomplete_colors" parameter sets
        whether showing these stars in a randomly generated color (True) or
        displaying them according to their stellar population if given (False).
        Caution: autocomplete_colors only avoid repetition if the colors are
        given by their hexadecimal name.
    independent : boolean
        Different figures for each UV and WV plane? Default is False.
    gf_uw2 : str
        The gf_uw2 parameter sets the aspect ratio of the Toomre diagram. The
        values, circumferences (default) or ellipses, indicate how lines of
        constant total velocity look in the resulting diagram.
    asso_par : str, optional
        Name (without extension) of the custom file for the parameters of the
        associations and moving groups. The CSV file must have the same
        structure as the default provided and be located in the current
        directory.

    """

    # Read the input file and convert to DataFrame
    if '.txt' in input_file:
        df = pd.read_csv(input_file, sep=' ', index_col=None)
    elif '.csv' in input_file:
        df = pd.read_csv(input_file, index_col=None)
    else:
        raise TypeError('Only txt or csv files are allowed.')

    # Check the names of the columns
    # The column names required are as follows:
    if 'Name' not in df.columns:
        raise Exception('Not found column called Name')
    if 'RA_J2000_deg' not in df.columns:
        raise Exception('Not found column called RA_J2000_deg')
    if 'DE_J2000_deg' not in df.columns:
        raise Exception('Not found column called DE_J2000_deg')
    if 'muRA_masa-1' not in df.columns:
        raise Exception('Not found column called muRA_masa-1')
    if 'emuRA_masa-1' not in df.columns:
        raise Exception('Not found column called emuRA_masa-1')
    if 'muDE_masa-1' not in df.columns:
        raise Exception('Not found column called muDE_masa-1')
    if 'emuDE_masa-1' not in df.columns:
        raise Exception('Not found column called emuDE_masa-1')
    if 'Vr_kms-1' not in df.columns:
        raise Exception('Not found column called Vr_kms-1')
    if 'eVr_kms-1' not in df.columns:
        raise Exception('Not found column called eVr_kms-1')
    if 'd_pc' not in df.columns and 'pi_mas' not in df.columns:
        raise Exception('Not found column called d_pc or pi_mas')
    elif ('d_pc' in df.columns and 'ed_pc' not in df.columns and
          'pi_mas' not in df.columns):
        raise Exception('Not found column called ed_pc')
    elif ('pi_mas' in df.columns and 'epi_mas' not in df.columns and
          'd_pc' not in df.columns):
        raise Exception('Not found column called epi_mas')

    # If PLX and/or EPLX is not given in the input file,
    # they are calculated from d and ed.
    if 'pi_mas' not in df.columns:
        print('Parallax and its error are calculated from d_pc and ed_pc')
        plx = []
        eplx = []
        for i in range(len(df['d_pc'])):
            plx.append(1000/df['d_pc'][i])
            eplx.append(1000*df['ed_pc'][i]/(df['d_pc'][i])**(2))
        # Add the calculated values to the DataFrame
        df['pi_mas'] = plx
        df['epi_mas'] = eplx
    else:
        pass

    # Handle zeros in the errors
    if errors == 'zeros':
        df['emuRA_masa-1'] = df['emuRA_masa-1'].fillna(0.)
        df['emuDE_masa-1'] = df['emuDE_masa-1'].fillna(0.)
        df['eVr_kms-1'] = df['eVr_kms-1'].fillna(0.)
        df['epi_mas'] = df['epi_mas'].fillna(0.)
    elif errors == 'drop':
        df['emuRA_masa-1'].replace(0., np.nan, inplace=True)
        df['emuDE_masa-1'].replace(0., np.nan, inplace=True)
        df['eVr_kms-1'].replace(0., np.nan, inplace=True)
        df['epi_mas'].replace(0., np.nan, inplace=True)
        df.dropna(subset=['epi_mas', 'emuRA_masa-1', 'emuDE_masa-1', 
                          'eVr_kms-1'], inplace=True)
        df.reset_index(drop=True, inplace=True)
    else:
        raise ValueError('Options are zeros or drop')

    df_out = spk_groups(df, file_name, asso_par)
    spk_graphs(df_out, file_name, dict_colors_groups,
               dict_colors_intergroup_stars, autocomplete_colors,
               independent, gf_uw2, asso_par)


class CommandLine:
    def __init__(self):
        parser = argparse.ArgumentParser(add_help=True)
        parser.add_argument("-i", help="the input file containing the data" +
                            " (.txt or .csv)", required=True, default="")
        parser.add_argument("-e", help="set how the code deals with zeros" +
                            " in the errors: drop (default) or zeros",
                            required=False, default="drop")
        parser.add_argument("-fn", help="set a string that all the output" +
                            " files must contain", required=False, default="")
        parser.add_argument("-dcg", help="set colors for new groups or" +
                            " change those given by default (str-like python" +
                            " dictionary)", required=False, default=None)
        parser.add_argument("-dci", help="set colors for new intergroup" +
                            " stars or change those given by default" +
                            " (str-like python dictionary)", required=False,
                            default=None)
        parser.add_argument("-ac", help="should the stars candidates to any" +
                            " association, group, or intergroup not defined" +
                            " in the dictionaries be shown in a random" +
                            " color? Default is False", type=ast.literal_eval,
                            required=False, default=False)
        parser.add_argument("-ind", help="Independent graphs for each plane?" +
                            " Default is False", type=ast.literal_eval,
                            required=False, default=False)
        parser.add_argument("-g", help="set the aspect ratio of the Toomre" +
                            " diagram: circumferences (default) or ellipses",
                            required=False, default="circumferences")
        parser.add_argument("-asso", help=" Name (without extension) of the" +
                            " custom file for the parameters of the" +
                            " associations and moving groups.",
                            required=False, default=None)
        argument = parser.parse_args()
        spk_wrapper(argument.i, argument.e, argument.fn, argument.dcg,
                    argument.dci, argument.ac, argument.ind, argument.g,
                    argument.asso)


if __name__ == '__main__':
    app = CommandLine()
