import os
import numpy as np
from math import atan
import pandas as pd


def spk_groups(df, file_name="", asso_par=None):

    def groups(ra_j2000, dec_j2000, mu_ra, e_mu_ra, mu_dec, e_mu_dec, pi, e_pi,
               rv, e_rv, asso_par=None):
        """
        Function to calculate the Galactic-space velocity components and their
        errors, evaluate membership to kinematic moving groups and
        associations, and assign stellar populations for each star.
        The Galactic-space velocity components and their errors are calculated
        following Johnson and Soderblom (1987), while the stellar populations
        are assigned according to Bensby et al. 2003, 2005.

        Parameters
        ----------

        ra_j2000, dec_j2000 : float
            Equatorial coordinates (degrees) of the star in J2000 and
            epoch 2000
        mu_ra, e_mu_ra, mu_dec, e_mu_dec: float
            Proper motions and their errors in mas yr-1
        pi, e_pi: float
            Parallax and its error in mas
        rv, e_rv: float
            Radial velocity and its error in km s-1
        asso_par : str, optional
            Name (without extension) of the custom file for the parameters of the
            associations and moving groups. The CSV file must have the same
            structure as the default provided and be located in the current
            directory.
        """

        # The equivalent in km s-1 of 1 AU in one tropical year.
        Av = 4.74047

        # Convert degrees to radians and mas to arcsec
        alfarad = ra_j2000*np.pi/180.
        deltarad = dec_j2000*np.pi/180.
        muacor = mu_ra/1000.
        emuacor = e_mu_ra/1000.
        mudcor = mu_dec/1000.
        emudcor = e_mu_dec/1000.
        picor = pi/1000.
        epicor = e_pi/1000.

        # Coordinates of the North Galatic Pole (ICRS(J2000))
        ra_ngp = 192.85948*np.pi/180
        dec_ngp = 27.12825*np.pi/180
        l_omega = 32.93192*np.pi/180

        # t matrix
        t = np.zeros((3, 3))

        t[0, 0] = (-np.sin(dec_ngp)*np.cos(ra_ngp)*np.cos((np.pi/2) +
                   l_omega)-np.sin(ra_ngp)*np.sin((np.pi/2)+l_omega))
        t[0, 1] = (-np.sin(dec_ngp)*np.sin(ra_ngp)*np.cos((np.pi/2)+l_omega) +
                   np.cos(ra_ngp)*np.sin((np.pi/2)+l_omega))
        t[0, 2] = np.cos(dec_ngp)*np.cos((np.pi/2)+l_omega)
        t[1, 0] = (-np.sin(dec_ngp)*np.cos(ra_ngp)*np.sin((np.pi/2)+l_omega) +
                   np.sin(ra_ngp)*np.cos((np.pi/2)+l_omega))
        t[1, 1] = (-np.sin(dec_ngp)*np.sin(ra_ngp)*np.sin((np.pi/2)+l_omega) -
                   np.cos(ra_ngp)*np.cos((np.pi/2)+l_omega))
        t[1, 2] = np.cos(dec_ngp)*np.sin((np.pi/2)+l_omega)
        t[2, 0] = np.cos(dec_ngp)*np.cos(ra_ngp)
        t[2, 1] = np.cos(dec_ngp)*np.sin(ra_ngp)
        t[2, 2] = np.sin(dec_ngp)

        # To optimize memory usage when working with large datasets, the user
        # can use the values commented, which are accurate enough for typical
        # errors in the parameters, instead of the expressions above.
        # t[0, 0] = -5.48755604*10**(-2)
        # t[0, 1] = -8.734370902*10**(-1)
        # t[0, 2] = -4.838350155*10**(-1)
        # t[1, 0] = 4.941094279*10**(-1)
        # t[1, 1] = -4.448296300*10**(-1)
        # t[1, 2] = 7.469822445*10**(-1)
        # t[2, 0] = -8.676661490*10**(-1)
        # t[2, 1] = -1.980763734*10**(-1)
        # t[2, 2] = 4.559837762*10**(-1)

        def M(alfa, delta):
            m = np.zeros((3, 3))
            m[0, 0] = np.cos(alfa)*np.cos(delta)
            m[0, 1] = -np.sin(alfa)
            m[0, 2] = -np.cos(alfa)*np.sin(delta)
            m[1, 0] = np.sin(alfa)*np.cos(delta)
            m[1, 1] = np.cos(alfa)
            m[1, 2] = -np.sin(alfa)*np.sin(delta)
            m[2, 0] = np.sin(delta)
            m[2, 1] = 0.
            m[2, 2] = np.cos(delta)
            return m

        def D(mu_alfa, mu_delta, rho, paralax, k):
            d1 = np.zeros(3)
            d1[0] = rho
            d1[1] = k*mu_alfa/paralax
            d1[2] = k*mu_delta/paralax
            return d1

        r = M(alfarad, deltarad)
        d = D(muacor, mudcor, rv, picor, Av)

        b = np.zeros((3, 3))
        for l in range(len(t)):
            for j in range(len(r[0])):
                for k in range(len(r)):
                    b[l, j] += t[l, k]*r[k, j]

        u = b[0, 0]*d[0]+b[0, 1]*d[1]+b[0, 2]*d[2]
        v = b[1, 0]*d[0]+b[1, 1]*d[1]+b[1, 2]*d[2]
        w = b[2, 0]*d[0]+b[2, 1]*d[1]+b[2, 2]*d[2]

        # Compute the errors assuming the variables are uncorrelated
        C = b**(2)
        s = np.zeros(3)
        s[0] = (e_rv)**(2)
        s[1] = (Av/picor)**(2)*(emuacor**(2)+(muacor*epicor/picor)**(2))
        s[2] = (Av/picor)**(2)*(emudcor**(2)+(mudcor*epicor/picor)**(2))

        le = np.zeros(3)
        le[0] = (2*muacor*mudcor*(Av)**(2)*(epicor)**(2) /
                 (picor)**(4))*b[0, 1]*b[0, 2]
        le[1] = (2*muacor*mudcor*(Av)**(2)*(epicor)**(2) /
                 (picor)**(4))*b[1, 1]*b[1, 2]
        le[2] = (2*muacor*mudcor*(Av)**(2)*(epicor)**(2) /
                 (picor)**(4))*b[2, 1]*b[2, 2]

        eucov = np.sqrt((C[0, 0]*s[0]+le[0])+(C[0, 1]*s[1]+le[1]) +
                        (C[0, 2]*s[2]+le[2]))
        evcov = np.sqrt((C[1, 0]*s[0]+le[0])+(C[1, 1]*s[1]+le[1]) +
                        (C[1, 2]*s[2]+le[2]))
        ewcov = np.sqrt((C[2, 0]*s[0]+le[0])+(C[2, 1]*s[1]+le[1]) +
                        (C[2, 2]*s[2]+le[2]))

        # Eggen's definition of the young disk (YD)
        u1 = -7.0 + ((v+30.0)/(-5.0+30.0))*(20.0+7.0)
        u2 = -28.0 + ((-48.5+28.0)/(-2.0+14.5))*(-2-v)
        u3 = 0.0 + ((v+2.0)/(2.7+2.0))*(8.0-0.0)

        # Does the star belong (YD) or not (NYD) to the young disk region as
        # defined by Eggen?
        if (-48.0 <= u <= 20.0 and -30.0 <= v <= -2.0 and u2 <= u <= u1 and
                -30.0 <= w <= 15.0):
            lab_YD = 'YD'
        elif u3 <= u <= 20.0 and -2.0 <= v <= 2.0 and -30.0 <= w <= 15.0:
            lab_YD = 'YD'
        else:
            lab_YD = 'NYD'

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

        # Ellipsoides in the Galactic-space velocity represent the associations
        # and moving groups
        def elipsoide(xg, yg, zg, aa, bb, cc, xs, ys, zs):
            dxy = np.sqrt((xs-xg)**(2)+(ys-yg)**(2))
            dt = np.sqrt((xs-xg)**(2) + (ys-yg)**(2) + (zs-zg)**(2))
            RA = atan((ys-yg)/(xs-xg))
            ZZ = abs(zs-zg)
            RZ = atan(ZZ/dxy)
            xe = aa*np.cos(RA)*np.cos(RZ)
            ye = bb*np.sin(RA)*np.cos(RZ)
            ze = cc*np.sin(RZ)
            dte = np.sqrt(xe**(2)+ye**(2)+ze**(2))
            return dt, dte

        dt = []
        dte = []
        for i in range(len(df_par)):
            dt_i, dte_i = elipsoide(df_par['U_kms-1'][i],
                                    df_par['V_kms-1'][i],
                                    df_par['W_kms-1'][i],
                                    df_par['sigU_kms-1'][i],
                                    df_par['sigV_kms-1'][i],
                                    df_par['sigW_kms-1'][i], u, v, w)
            dt.append(dt_i)
            dte.append(dte_i)

        # Is the star candidate to any of the moving groups or associations
        # defined in the parameter file association_parameters.csv?
        lab_skg = None
        for i in range(len(dt)):
            if dt[i] <= dte[i]:
                if lab_skg is None:
                    lab_skg = df_par['ASSOCIATION/GROUP'][i]
                else:
                    lab_skg += '/' + df_par['ASSOCIATION/GROUP'][i]
            if i == len(df_par)-1 and lab_skg is None:
                lab_skg = 'NNN'

        # Correction to a dynamical local standard of rest (LSR):
        # The Sun's peculiar motion relative to the LSR is assumed to be
        # (U, V, W) = (+10.00,+5.25,+7.17) km s-1 (Dehnen & Binney 1998)
        # A comprehensive list of values of the solar peculiar velocity
        # components can be found in Francis and Anderson 2009.
        u_s = +10.00
        v_s = +5.25
        w_s = +7.17
        u_LSR = u + u_s
        v_LSR = v + v_s
        w_LSR = w + w_s

        def prob(u_LSR, sig_u, v_LSR, v_asym, sig_v, w_LSR, sig_w, X_ns):
            """
            Function to calculate the probability that a star belongs to a
            certain stellar population

            Parameters
            ----------

            u_LSR, v_LSR, w_LSR : float
                The Galactic-space velocity components corrected for the
                peculiar motion of the Sun relative to the LSR.
            sig_u, sig_v, sig_w : float
                The characteristic velocity dispersions of the population
            v_asym : float
                The asymmetric drift
            X_ns : float
                The assumed fraction of stars of each population in the
                solar neighborhood
            """

            k_nor = 1.0/(((2.0*np.pi)**(3.0/2.0))*sig_u*sig_v*sig_w)
            u_p = (u_LSR**(2.0))/(2.0*(sig_u**(2.0)))
            v_p = ((v_LSR-v_asym)**(2.0))/(2.0*(sig_v**(2.0)))
            w_p = (w_LSR**(2.0))/(2.0*(sig_w**(2.0)))
            prob1 = X_ns*k_nor*(np.exp(-u_p-v_p-w_p))
            return prob1

        # CSV file containing the parameters to calculate the probabilities
        try:
            df_par2 = pd.read_csv(
                os.path.join(os.path.dirname(__file__),
                             'param_prob_populations.csv'), index_col=None)
        except FileNotFoundError as e:
            raise FileNotFoundError('The file param_prob_populations.csv ' +
                                    'was not found in the corresponding ' +
                                    'directory.').with_traceback(e.__traceback__)

        # Thin disk (D) probabilities
        p_D = prob(u_LSR, df_par2['sig_u_kms-1'][0],
                   v_LSR, df_par2['v_asym_kms-1'][0],
                   df_par2['sig_v_kms-1'][0], w_LSR,
                   df_par2['sig_w_kms-1'][0],
                   df_par2['X_ns'][0])

        # Thick disk (TD) probabilities
        p_TD = prob(u_LSR, df_par2['sig_u_kms-1'][1],
                    v_LSR, df_par2['v_asym_kms-1'][1],
                    df_par2['sig_v_kms-1'][1],
                    w_LSR, df_par2['sig_w_kms-1'][1],
                    df_par2['X_ns'][1])

        # Halo (H) probabilities
        p_H = prob(u_LSR, df_par2['sig_u_kms-1'][2],
                   v_LSR, df_par2['v_asym_kms-1'][2],
                   df_par2['sig_v_kms-1'][2], w_LSR,
                   df_par2['sig_w_kms-1'][2],
                   df_par2['X_ns'][2])

        # Relative probabilities
        p_TD_D = p_TD/p_D
        p_TD_H = p_TD/p_H

        #  Population assignment
        if p_TD_H > 1.0:
            if p_TD_D <= 0.6:
                lab_kin = 'D'
            elif p_TD_D > 0.6 and p_TD_D < 2.0:
                lab_kin = 'TD-D'
            elif p_TD_D >= 2.0:
                lab_kin = 'TD'
        else:
            lab_kin = 'H'

        return lab_YD, lab_skg, lab_kin, u, eucov, v, evcov, w, ewcov

    lab_YD, lab_skg, lab_kin, u, eucov, v, evcov, w, ewcov = (
        [] for i in range(9))

    for i in range(len(df['RA_J2000_deg'])):
        yd, skg, kin, u_g, eu_g, v_g, ev_g, w_g, ew_g = (
            groups(df['RA_J2000_deg'][i], df['DE_J2000_deg'][i],
                   df['muRA_masa-1'][i], df['emuRA_masa-1'][i],
                   df['muDE_masa-1'][i], df['emuDE_masa-1'][i],
                   df['pi_mas'][i], df['epi_mas'][i],
                   df['Vr_kms-1'][i], df['eVr_kms-1'][i]))
        lab_YD.append(yd)
        lab_skg.append(skg)
        lab_kin.append(kin)
        u.append(u_g)
        eucov.append(eu_g)
        v.append(v_g)
        evcov.append(ev_g)
        w.append(w_g)
        ewcov.append(ew_g)

    dict1 = {'Name': df['Name'],
              'Vr_kms-1': df['Vr_kms-1'], 'eVr_kms-1': df['eVr_kms-1'],
              'U_kms-1': u, 'eU_kms-1': eucov, 'V_kms-1': v, 'eV_kms-1': evcov,
              'W_kms-1': w, 'eW_kms-1': ewcov, 'YD': lab_YD,
              'Kin_Group': lab_skg, 'Kin_Pop': lab_kin}
    df_out = pd.DataFrame(dict1, columns=dict1.keys())

    formats = {'Vr_kms-1': "{:.5f}", 'eVr_kms-1': "{:.5f}", 
                'U_kms-1': "{:.5f}", 'eU_kms-1': "{:.5f}",
                'V_kms-1': "{:.5f}", 'eV_kms-1': "{:.5f}",
                'W_kms-1': "{:.5f}", 'eW_kms-1': "{:.5f}"}

    for col, f in formats.items():
        df_out[col] = df_out[col].map(lambda x: f.format(x))

    df_out.set_index('Name')

    out_name = 'SteParKin_'+file_name+'.csv'
    df_out.to_csv(out_name, sep=',', index=False)

    return df_out
