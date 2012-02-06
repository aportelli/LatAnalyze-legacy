/* latan_nunits.h, part of LatAnalyze library
 *
 * Copyright (C) 2010, 2011, 2012 Antonin Portelli
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef LATAN_NUNITS_H_
#define LATAN_NUNITS_H_

/* sources : Review of Particle Physics 2010 (PDG)  ( http://pdglive.lbl.gov/ )
 *           Review of Lattice Results  2010 (FLAG) ( arXiv:1011.4408v1 )
 */

/* units */
#define NU_FM               0.005067731             /* MeV^(-1) */

/* constants */
#define NU_ALPHA_EM         7.2973525376e-3         /* no dim */
#define NU_ELECTRON_CHARGE  0.302822120214353       /* no dim */

/* masses */
/** pi **/
#define NU_M_pi_p           139.57018               /* MeV */
#define NU_M_pi_p_ERR       3.5e-4                  /* MeV */
#define NU_M_pi_0           134.9766                /* MeV */
#define NU_M_pi_0_ERR       6.0e-4                  /* MeV */
#define NU_M_pi_m           139.57018               /* MeV */
#define NU_M_pi_m_ERR       3.5e-4                  /* MeV */
#define NU_M_pi_iso         134.8                   /* MeV */
#define NU_M_pi_iso_ERR     0.3                     /* MeV */
#define NU_M_pi             NU_M_pi_iso
#define NU_M_pi_ERR         NU_M_pi_iso_ERR
#define NU_M_pi_p_miso      139.5                   /* MeV */
#define NU_M_pi_p_miso_ERR  0.3                     /* MeV */
#define NU_M_pi_0_miso      135.1                   /* MeV */
#define NU_M_pi_0_miso_ERR  0.3                     /* MeV */

/** K **/
#define NU_M_K_p            493.677                 /* MeV */
#define NU_M_K_p_ERR        0.016                   /* MeV */
#define NU_M_K_0            497.614                 /* MeV */
#define NU_M_K_0_ERR        0.024                   /* MeV */
#define NU_M_K_m            493.677                 /* MeV */
#define NU_M_K_m_ERR        0.016                   /* MeV */
#define NU_M_K_iso          494.2                   /* MeV */
#define NU_M_K_iso_ERR      0.5                     /* MeV */
#define NU_M_Kchi           486.3733                /* MeV */
#define NU_M_Kchi_ERR       0.01                    /* MeV */ /* TODO : real error */

/** rho **/
#define NU_M_rho_p          775.49                  /* MeV */
#define NU_M_rho_p_ERR      0.34                    /* MeV */
#define NU_M_rho_0          775.49                  /* MeV */
#define NU_M_rho_0_ERR      0.34                    /* MeV */
#define NU_M_rho_m          775.49                  /* MeV */
#define NU_M_rho_m_ERR      0.34                    /* MeV */
#define NU_M_rho            NU_M_rho_0
#define NU_M_rho_ERR        NU_M_rho_0_ERR

/** K* (892) **/
#define NU_M_Kst_p          891.66                  /* MeV */
#define NU_M_Kst_p_ERR      0.26                    /* MeV */
#define NU_M_Kst_0          895.94                  /* MeV */
#define NU_M_Kst_0_ERR      0.22                    /* MeV */
#define NU_M_Kst_m          891.66                  /* MeV */
#define NU_M_Kst_m_ERR      0.26                    /* MeV */
#define NU_M_Kst            893.8                   /* MeV */
#define NU_M_Kst_ERR        0.17                    /* MeV */

/** N **/
#define NU_M_p              938.272013              /* MeV */
#define NU_M_p_ERR          0.000023                /* MeV */
#define NU_M_n              939.565346              /* MeV */
#define NU_M_n_ERR          0.000023                /* MeV */
#define NU_M_N              938.9186795             /* MeV */
#define NU_M_N_ERR          0.0000115               /* MeV */

/** Sigma **/
#define NU_M_Sigma_p        1189.37                 /* MeV */
#define NU_M_Sigma_p_ERR    0.07                    /* MeV */
#define NU_M_Sigma_0        1192.642                /* MeV */
#define NU_M_Sigma_0_ERR    0.024                   /* MeV */
#define NU_M_Sigma_m        1197.449                /* MeV */
#define NU_M_Sigma_m_ERR    0.03                    /* MeV */
#define NU_M_Sigma          1193.154                /* MeV */
#define NU_M_Sigma_ERR      0.027                   /* MeV */

/** Xi **/
#define NU_M_Xi_0           1314.86                 /* MeV */
#define NU_M_Xi_0_ERR       0.2                     /* MeV */
#define NU_M_Xi_m           1321.71                 /* MeV */
#define NU_M_Xi_m_ERR       0.07                    /* MeV */
#define NU_M_Xi             1318.285                /* MeV */
#define NU_M_Xi_ERR         0.11                    /* MeV */

/** Delta (1232) **/
#define NU_M_Delta_pp       1232.0                  /* MeV */
#define NU_M_Delta_pp_ERR   1.0                     /* MeV */
#define NU_M_Delta_p        1232.0                  /* MeV */
#define NU_M_Delta_p_ERR    1.0                     /* MeV */
#define NU_M_Delta_0        1232.0                  /* MeV */
#define NU_M_Delta_0_ERR    1.0                     /* MeV */
#define NU_M_Delta_m        1232.0                  /* MeV */
#define NU_M_Delta_m_ERR    1.0                     /* MeV */
#define NU_M_Delta          NU_M_Delta_0
#define NU_M_Delta_ERR      NU_M_Delta_0_ERR

/** Sigma* (1385) **/
#define NU_M_Sigmast_p      1382.8                  /* MeV */
#define NU_M_Sigmast_p_ERR  0.4                     /* MeV */
#define NU_M_Sigmast_0      1383.7                  /* MeV */
#define NU_M_Sigmast_0_ERR  1.0                     /* MeV */
#define NU_M_Sigmast_m      1387.2                  /* MeV */
#define NU_M_Sigmast_m_ERR  0.5                     /* MeV */
#define NU_M_Sigmast        1384.57                 /* MeV */
#define NU_M_Sigmast_ERR    0.40                    /* MeV */

/** Xi* (1530) **/
#define NU_M_Xist_0         1531.8                  /* MeV */
#define NU_M_Xist_0_ERR     0.32                    /* MeV */
#define NU_M_Xist_m         1535.0                  /* MeV */
#define NU_M_Xist_m_ERR     0.6                     /* MeV */
#define NU_M_Xist           1533.4                  /* MeV */
#define NU_M_Xist_ERR       0.34                    /* MeV */

/** Omega **/
#define NU_M_Omega_m        1672.45                 /* MeV */
#define NU_M_Omega_m_ERR    0.29                    /* MeV */
#define NU_M_Omega          NU_M_Omega_m
#define NU_M_Omega_ERR      NU_M_Omega_m_ERR

#endif
