/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *--------------------------------------------------------------*/

#include <arkode/arkode.h>
#include <arkode/arkode_sprk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode/arkode_butcher.h"
#include "arkode_impl.h"

static ARKodeSPRKTable arkodeSymplecticEuler(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(1);
  if (!sprk_table) { return NULL; }
  sprk_table->q       = 1;
  sprk_table->stages  = 1;
  sprk_table->a[0]    = SUN_RCONST(1.0);
  sprk_table->ahat[0] = SUN_RCONST(1.0);
  return sprk_table;
}

/*
  The following methods are from:

  J Candy, W Rozmus, A symplectic integration algorithm for separable
  Hamiltonian functions, Journal of Computational Physics, Volume 92, Issue 1,
  1991, Pages 230-256, ISSN 0021-9991,
  https://doi.org/10.1016/0021-9991(91)90299-Z.
 */

static ARKodeSPRKTable arkodeSymplecticLeapfrog2(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(2);
  if (!sprk_table) { return NULL; }
  sprk_table->q       = 2;
  sprk_table->stages  = 2;
  sprk_table->a[0]    = SUN_RCONST(0.5);
  sprk_table->a[1]    = SUN_RCONST(0.5);
  sprk_table->ahat[0] = SUN_RCONST(0.0);
  sprk_table->ahat[1] = SUN_RCONST(1.0);
  return sprk_table;
}

static ARKodeSPRKTable arkodeSymplecticPseudoLeapfrog2(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(2);
  if (!sprk_table) { return NULL; }
  sprk_table->q       = 2;
  sprk_table->stages  = 2;
  sprk_table->a[0]    = SUN_RCONST(1.0);
  sprk_table->a[1]    = SUN_RCONST(0.0);
  sprk_table->ahat[0] = SUN_RCONST(0.5);
  sprk_table->ahat[1] = SUN_RCONST(0.5);
  return sprk_table;
}

static ARKodeSPRKTable arkodeSymplecticCandyRozmus4(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(4);
  if (!sprk_table) { return NULL; }
  sprk_table->q      = 4;
  sprk_table->stages = 4;
  sprk_table->a[0] =
    (SUN_RCONST(2.0) +
     SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0) / SUN_RCONST(3.0)) +
     SUNRpowerR(SUN_RCONST(2.0), -SUN_RCONST(1.0) / SUN_RCONST(3.0))) /
    SUN_RCONST(6.0);
  sprk_table->a[1] =
    (SUN_RCONST(1.0) -
     SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0) / SUN_RCONST(3.0)) -
     SUNRpowerR(SUN_RCONST(2.0), -SUN_RCONST(1.0) / SUN_RCONST(3.0))) /
    SUN_RCONST(6.0);
  sprk_table->a[2]    = sprk_table->a[1];
  sprk_table->a[3]    = sprk_table->a[0];
  sprk_table->ahat[0] = SUN_RCONST(0.0);
  sprk_table->ahat[1] =
    SUN_RCONST(1.0) /
    (SUN_RCONST(2.0) -
     SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0) / SUN_RCONST(3.0)));
  sprk_table->ahat[2] =
    SUN_RCONST(1.0) /
    (SUN_RCONST(1.0) -
     SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(2.0) / SUN_RCONST(3.0)));
  sprk_table->ahat[3] = sprk_table->ahat[1];
  return sprk_table;
}

/*
  The following methods are from:

  Ruth, R. D. (1983). A CANONICAL INTEGRATION TECHNIQUE.
  IEEE Transactions on Nuclear Science, 30(4).
  https://accelconf.web.cern.ch/p83/PDF/PAC1983_2669.PDF
 */

static ARKodeSPRKTable arkodeSymplecticRuth3(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(3);
  if (!sprk_table) { return NULL; }
  sprk_table->q       = 3;
  sprk_table->stages  = 3;
  sprk_table->a[0]    = SUN_RCONST(2.0) / SUN_RCONST(3.0);
  sprk_table->a[1]    = -SUN_RCONST(2.0) / SUN_RCONST(3.0);
  sprk_table->a[2]    = SUN_RCONST(1.0);
  sprk_table->ahat[0] = SUN_RCONST(7.0) / SUN_RCONST(24.0);
  sprk_table->ahat[1] = SUN_RCONST(3.0) / SUN_RCONST(4.0);
  sprk_table->ahat[2] = -SUN_RCONST(1.0) / SUN_RCONST(24.0);
  return sprk_table;
}

/*
  The following methods are from:

  McLachlan, R.I., Atela, P.: The accuracy of symplectic integrators.
  Nonlinearity. 5, 541–562 (1992). https://doi.org/10.1088/0951-7715/5/2/011
 */

static ARKodeSPRKTable arkodeSymplecticMcLachlan2(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(2);
  if (!sprk_table) { return NULL; }
  sprk_table->q      = 2;
  sprk_table->stages = 2;
  sprk_table->a[1]   = SUN_RCONST(1.0) - (SUN_RCONST(1.0) / SUN_RCONST(2.0)) *
                                         SUNRsqrt(SUN_RCONST(2.0));
  sprk_table->a[0]    = SUN_RCONST(1.0) - sprk_table->a[1];
  sprk_table->ahat[1] = SUN_RCONST(1.0) /
                        (SUN_RCONST(2.0) * (SUN_RCONST(1.0) - sprk_table->a[1]));
  sprk_table->ahat[0] = SUN_RCONST(1.0) - sprk_table->ahat[1];
  return sprk_table;
}

static ARKodeSPRKTable arkodeSymplecticMcLachlan3(void)
{
  sunrealtype w              = SUN_RCONST(0.0);
  sunrealtype y              = SUN_RCONST(0.0);
  sunrealtype z              = SUN_RCONST(0.0);
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(3);
  if (!sprk_table) { return NULL; }

  sprk_table->q      = 3;
  sprk_table->stages = 3;

  z = -SUNRpowerR((SUN_RCONST(2.0) / SUN_RCONST(27.0)) -
                    SUN_RCONST(1.0) /
                      (SUN_RCONST(9.0) * SUNRsqrt(SUN_RCONST(3.0))),
                  SUN_RCONST(1.0) / SUN_RCONST(3.0));
  w = -SUN_RCONST(2.0) / SUN_RCONST(3.0) +
      SUN_RCONST(1.0) / (SUN_RCONST(9.0) * z) + z;
  y                = (SUN_RCONST(1.0) + w * w) / SUN_RCONST(4.0);
  sprk_table->a[0] = SUNRsqrt(SUN_RCONST(1.0) / (SUN_RCONST(9.0) * y) -
                              w / SUN_RCONST(2.0) + SUNRsqrt(y)) -
                     SUN_RCONST(1.0) / (SUN_RCONST(3.0) * SUNRsqrt(y));
  sprk_table->a[1] = SUN_RCONST(0.25) / sprk_table->a[0] -
                     sprk_table->a[0] / SUN_RCONST(2.0);
  sprk_table->a[2]    = SUN_RCONST(1.0) - sprk_table->a[0] - sprk_table->a[1];
  sprk_table->ahat[0] = sprk_table->a[2];
  sprk_table->ahat[1] = sprk_table->a[1];
  sprk_table->ahat[2] = sprk_table->a[0];
  return sprk_table;
}

static ARKodeSPRKTable arkodeSymplecticMcLachlan4(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(4);
  if (!sprk_table) { return NULL; }
  sprk_table->q       = 4;
  sprk_table->stages  = 4;
  sprk_table->a[0]    = SUN_RCONST(0.515352837431122936);
  sprk_table->a[1]    = -SUN_RCONST(0.085782019412973646);
  sprk_table->a[2]    = SUN_RCONST(0.441583023616466524);
  sprk_table->a[3]    = SUN_RCONST(0.128846158365384185);
  sprk_table->ahat[0] = SUN_RCONST(0.134496199277431089);
  sprk_table->ahat[1] = -SUN_RCONST(0.224819803079420806);
  sprk_table->ahat[2] = SUN_RCONST(0.756320000515668291);
  sprk_table->ahat[3] = SUN_RCONST(0.33400360328632142);
  return sprk_table;
}

static ARKodeSPRKTable arkodeSymplecticMcLachlan5(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(6);
  if (!sprk_table) { return NULL; }
  sprk_table->q       = 5;
  sprk_table->stages  = 6;
  sprk_table->a[0]    = SUN_RCONST(0.339839625839110000);
  sprk_table->a[1]    = -SUN_RCONST(0.088601336903027329);
  sprk_table->a[2]    = SUN_RCONST(0.5858564768259621188);
  sprk_table->a[3]    = -SUN_RCONST(0.603039356536491888);
  sprk_table->a[4]    = SUN_RCONST(0.3235807965546976394);
  sprk_table->a[5]    = SUN_RCONST(0.4423637942197494587);
  sprk_table->ahat[0] = SUN_RCONST(0.1193900292875672758);
  sprk_table->ahat[1] = SUN_RCONST(0.6989273703824752308);
  sprk_table->ahat[2] = -SUN_RCONST(0.1713123582716007754);
  sprk_table->ahat[3] = SUN_RCONST(0.4012695022513534480);
  sprk_table->ahat[4] = SUN_RCONST(0.0107050818482359840);
  sprk_table->ahat[5] = -SUN_RCONST(0.0589796254980311632);
  return sprk_table;
}

/*
  The following methods are from:

  Yoshida, H.: Construction of higher order symplectic integrators.
  Phys Lett A. 150, 262–268 (1990).
  https://doi.org/10.1016/0375-9601(90)90092-3

 */

static ARKodeSPRKTable arkodeSymplecticYoshida6(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(8);
  if (!sprk_table) { return NULL; }
  sprk_table->q       = 6;
  sprk_table->stages  = 8;
  sprk_table->a[0]    = SUN_RCONST(0.7845136104775572638194976338663498757768);
  sprk_table->a[1]    = SUN_RCONST(0.2355732133593581336847931829785346016865);
  sprk_table->a[2]    = -SUN_RCONST(1.177679984178871006946415680964315734639);
  sprk_table->a[3]    = SUN_RCONST(1.315186320683911218884249728238862514352);
  sprk_table->a[4]    = sprk_table->a[2];
  sprk_table->a[5]    = sprk_table->a[1];
  sprk_table->a[6]    = sprk_table->a[0];
  sprk_table->a[7]    = SUN_RCONST(0.0);
  sprk_table->ahat[0] = sprk_table->a[0] / SUN_RCONST(2.0);
  sprk_table->ahat[1] = (sprk_table->a[0] + sprk_table->a[1]) / SUN_RCONST(2.0);
  sprk_table->ahat[2] = (sprk_table->a[1] + sprk_table->a[2]) / SUN_RCONST(2.0);
  sprk_table->ahat[3] = (sprk_table->a[2] + sprk_table->a[3]) / SUN_RCONST(2.0);
  sprk_table->ahat[4] = sprk_table->ahat[3];
  sprk_table->ahat[5] = sprk_table->ahat[2];
  sprk_table->ahat[6] = sprk_table->ahat[1];
  sprk_table->ahat[7] = sprk_table->ahat[0];
  return sprk_table;
}

/*
  The following methods are from:

  (Original) Suzuki, M., & Umeno, K. (1993). Higher-order decomposition theory
  of exponential operators and its applications to QMC and nonlinear dynamics.
  Computer simulation studies in condensed-matter physics VI, 74-86.
  https://doi.org/10.1007/978-3-642-78448-4_7

  McLachlan, R.I.: On the Numerical Integration of Ordinary Differential
  Equations by Symmetric Composition Methods. Siam J Sci Comput. 16, 151–168
  (1995). https://doi.org/10.1137/0916010

 */

static ARKodeSPRKTable arkodeSymplecticSuzukiUmeno816(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(16);
  if (!sprk_table) { return NULL; }
  sprk_table->q       = 8;
  sprk_table->stages  = 16;
  sprk_table->a[0]    = SUN_RCONST(0.7416703643506129534482278017838063156035);
  sprk_table->a[1]    = -SUN_RCONST(0.4091008258000315939973000958935634173099);
  sprk_table->a[2]    = SUN_RCONST(0.1907547102962383799538762564503716627355);
  sprk_table->a[3]    = -SUN_RCONST(0.5738624711160822666563877266355357421595);
  sprk_table->a[4]    = SUN_RCONST(0.2990641813036559238444635406886029882258);
  sprk_table->a[5]    = SUN_RCONST(0.3346249182452981837849579798821822886337);
  sprk_table->a[6]    = SUN_RCONST(0.3152930923967665966320566638110024309941);
  sprk_table->a[7]    = -SUN_RCONST(0.7968879393529163540197888401737330534463);
  sprk_table->a[8]    = sprk_table->a[6];
  sprk_table->a[9]    = sprk_table->a[5];
  sprk_table->a[10]   = sprk_table->a[4];
  sprk_table->a[11]   = sprk_table->a[3];
  sprk_table->a[12]   = sprk_table->a[2];
  sprk_table->a[13]   = sprk_table->a[1];
  sprk_table->a[14]   = sprk_table->a[0];
  sprk_table->a[15]   = SUN_RCONST(0.0);
  sprk_table->ahat[0] = sprk_table->a[0] / SUN_RCONST(2.0);
  sprk_table->ahat[1] = (sprk_table->a[0] + sprk_table->a[1]) / SUN_RCONST(2.0);
  sprk_table->ahat[2] = (sprk_table->a[1] + sprk_table->a[2]) / SUN_RCONST(2.0);
  sprk_table->ahat[3] = (sprk_table->a[2] + sprk_table->a[3]) / SUN_RCONST(2.0);
  sprk_table->ahat[4] = (sprk_table->a[3] + sprk_table->a[4]) / SUN_RCONST(2.0);
  sprk_table->ahat[5] = (sprk_table->a[4] + sprk_table->a[5]) / SUN_RCONST(2.0);
  sprk_table->ahat[6] = (sprk_table->a[5] + sprk_table->a[6]) / SUN_RCONST(2.0);
  sprk_table->ahat[7] = (sprk_table->a[6] + sprk_table->a[7]) / SUN_RCONST(2.0);
  sprk_table->ahat[8] = sprk_table->ahat[7];
  sprk_table->ahat[9] = sprk_table->ahat[6];
  sprk_table->ahat[10] = sprk_table->ahat[5];
  sprk_table->ahat[11] = sprk_table->ahat[4];
  sprk_table->ahat[12] = sprk_table->ahat[3];
  sprk_table->ahat[13] = sprk_table->ahat[2];
  sprk_table->ahat[14] = sprk_table->ahat[1];
  sprk_table->ahat[15] = sprk_table->ahat[0];
  return sprk_table;
}

/*
  The following methods are from:

  Sofroniou, M., Spaletta, G.: Derivation of symmetric composition constants for
  symmetric integrators. Optim Methods Softw. 20, 597–613 (2005).
  https://doi.org/10.1080/10556780500140664

 */

static ARKodeSPRKTable arkodeSymplecticSofroniou10(void)
{
  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(36);
  if (!sprk_table) { return NULL; }
  sprk_table->q      = 10;
  sprk_table->stages = 36;

  sprk_table->a[0]    = SUN_RCONST(0.078795722521686419263907679337684);
  sprk_table->a[1]    = SUN_RCONST(0.31309610341510852776481247192647);
  sprk_table->a[2]    = SUN_RCONST(0.027918383235078066109520273275299);
  sprk_table->a[3]    = -SUN_RCONST(0.22959284159390709415121339679655);
  sprk_table->a[4]    = SUN_RCONST(0.13096206107716486317465685927961);
  sprk_table->a[5]    = -SUN_RCONST(0.26973340565451071434460973222411);
  sprk_table->a[6]    = SUN_RCONST(0.074973343155891435666137105641410);
  sprk_table->a[7]    = SUN_RCONST(0.11199342399981020488957508073640);
  sprk_table->a[8]    = SUN_RCONST(0.36613344954622675119314812353150);
  sprk_table->a[9]    = -SUN_RCONST(0.39910563013603589787862981058340);
  sprk_table->a[10]   = SUN_RCONST(0.10308739852747107731580277001372);
  sprk_table->a[11]   = SUN_RCONST(0.41143087395589023782070411897608);
  sprk_table->a[12]   = -SUN_RCONST(0.0048663605831352617621956593099771);
  sprk_table->a[13]   = -SUN_RCONST(0.39203335370863990644808193642610);
  sprk_table->a[14]   = SUN_RCONST(0.051942502962449647037182904015976);
  sprk_table->a[15]   = SUN_RCONST(0.050665090759924496335874344156866);
  sprk_table->a[16]   = SUN_RCONST(0.049674370639729879054568800279461);
  sprk_table->a[17]   = SUN_RCONST(0.049317735759594537917680008339338);
  sprk_table->a[18]   = sprk_table->a[16];
  sprk_table->a[19]   = sprk_table->a[15];
  sprk_table->a[20]   = sprk_table->a[14];
  sprk_table->a[21]   = sprk_table->a[13];
  sprk_table->a[22]   = sprk_table->a[12];
  sprk_table->a[23]   = sprk_table->a[11];
  sprk_table->a[24]   = sprk_table->a[10];
  sprk_table->a[25]   = sprk_table->a[9];
  sprk_table->a[26]   = sprk_table->a[8];
  sprk_table->a[27]   = sprk_table->a[7];
  sprk_table->a[28]   = sprk_table->a[6];
  sprk_table->a[29]   = sprk_table->a[5];
  sprk_table->a[30]   = sprk_table->a[4];
  sprk_table->a[31]   = sprk_table->a[3];
  sprk_table->a[32]   = sprk_table->a[2];
  sprk_table->a[33]   = sprk_table->a[1];
  sprk_table->a[34]   = sprk_table->a[0];
  sprk_table->a[35]   = SUN_RCONST(0.0);
  sprk_table->ahat[0] = sprk_table->a[0] / SUN_RCONST(2.0);
  sprk_table->ahat[1] = (sprk_table->a[0] + sprk_table->a[1]) / SUN_RCONST(2.0);
  sprk_table->ahat[2] = (sprk_table->a[1] + sprk_table->a[2]) / SUN_RCONST(2.0);
  sprk_table->ahat[3] = (sprk_table->a[2] + sprk_table->a[3]) / SUN_RCONST(2.0);
  sprk_table->ahat[4] = (sprk_table->a[3] + sprk_table->a[4]) / SUN_RCONST(2.0);
  sprk_table->ahat[5] = (sprk_table->a[4] + sprk_table->a[5]) / SUN_RCONST(2.0);
  sprk_table->ahat[6] = (sprk_table->a[5] + sprk_table->a[6]) / SUN_RCONST(2.0);
  sprk_table->ahat[7] = (sprk_table->a[6] + sprk_table->a[7]) / SUN_RCONST(2.0);
  sprk_table->ahat[8] = (sprk_table->a[7] + sprk_table->a[8]) / SUN_RCONST(2.0);
  sprk_table->ahat[9] = (sprk_table->a[8] + sprk_table->a[9]) / SUN_RCONST(2.0);
  sprk_table->ahat[10] = (sprk_table->a[9] + sprk_table->a[10]) / SUN_RCONST(2.0);
  sprk_table->ahat[11] = (sprk_table->a[10] + sprk_table->a[11]) /
                         SUN_RCONST(2.0);
  sprk_table->ahat[12] = (sprk_table->a[11] + sprk_table->a[12]) /
                         SUN_RCONST(2.0);
  sprk_table->ahat[13] = (sprk_table->a[12] + sprk_table->a[13]) /
                         SUN_RCONST(2.0);
  sprk_table->ahat[14] = (sprk_table->a[13] + sprk_table->a[14]) /
                         SUN_RCONST(2.0);
  sprk_table->ahat[15] = (sprk_table->a[14] + sprk_table->a[15]) /
                         SUN_RCONST(2.0);
  sprk_table->ahat[16] = (sprk_table->a[15] + sprk_table->a[16]) /
                         SUN_RCONST(2.0);
  sprk_table->ahat[17] = (sprk_table->a[16] + sprk_table->a[17]) /
                         SUN_RCONST(2.0);
  sprk_table->ahat[18] = sprk_table->ahat[17];
  sprk_table->ahat[19] = sprk_table->ahat[16];
  sprk_table->ahat[20] = sprk_table->ahat[15];
  sprk_table->ahat[21] = sprk_table->ahat[14];
  sprk_table->ahat[22] = sprk_table->ahat[13];
  sprk_table->ahat[23] = sprk_table->ahat[12];
  sprk_table->ahat[24] = sprk_table->ahat[11];
  sprk_table->ahat[25] = sprk_table->ahat[10];
  sprk_table->ahat[26] = sprk_table->ahat[9];
  sprk_table->ahat[27] = sprk_table->ahat[8];
  sprk_table->ahat[28] = sprk_table->ahat[7];
  sprk_table->ahat[29] = sprk_table->ahat[6];
  sprk_table->ahat[30] = sprk_table->ahat[5];
  sprk_table->ahat[31] = sprk_table->ahat[4];
  sprk_table->ahat[32] = sprk_table->ahat[3];
  sprk_table->ahat[33] = sprk_table->ahat[2];
  sprk_table->ahat[34] = sprk_table->ahat[1];
  sprk_table->ahat[35] = sprk_table->ahat[0];

  return sprk_table;
}

ARKodeSPRKTable ARKodeSPRKTable_Create(int s, int q, const sunrealtype* a,
                                       const sunrealtype* ahat)
{
  if (s < 1 || !a || !ahat) { return NULL; }

  ARKodeSPRKTable sprk_table = ARKodeSPRKTable_Alloc(s);
  if (!sprk_table) { return NULL; }

  sprk_table->stages = s;
  sprk_table->q      = q;

  for (int i = 0; i < s; i++)
  {
    sprk_table->a[i]    = a[i];
    sprk_table->ahat[i] = ahat[i];
  }

  return sprk_table;
}

ARKodeSPRKTable ARKodeSPRKTable_Alloc(int stages)
{
  ARKodeSPRKTable sprk_table = NULL;

  sprk_table = (ARKodeSPRKTable)malloc(sizeof(struct ARKodeSPRKTableMem));
  if (!sprk_table) { return NULL; }

  memset(sprk_table, 0, sizeof(struct ARKodeSPRKTableMem));

  sprk_table->ahat = (sunrealtype*)malloc(stages * sizeof(sunrealtype));
  if (!(sprk_table->ahat))
  {
    ARKodeSPRKTable_Free(sprk_table);
    return NULL;
  }

  sprk_table->a = (sunrealtype*)malloc(stages * sizeof(sunrealtype));
  if (!(sprk_table->a))
  {
    ARKodeSPRKTable_Free(sprk_table);
    return NULL;
  }

  sprk_table->stages = stages;

  return sprk_table;
}

ARKodeSPRKTable ARKodeSPRKTable_Load(ARKODE_SPRKMethodID id)
{
  switch (id)
  {
  case ARKODE_SPRK_EULER_1_1: return arkodeSymplecticEuler();
  case ARKODE_SPRK_LEAPFROG_2_2: return arkodeSymplecticLeapfrog2();
  case ARKODE_SPRK_PSEUDO_LEAPFROG_2_2:
    return arkodeSymplecticPseudoLeapfrog2();
  case ARKODE_SPRK_RUTH_3_3: return arkodeSymplecticRuth3();
  case ARKODE_SPRK_MCLACHLAN_2_2: return arkodeSymplecticMcLachlan2();
  case ARKODE_SPRK_MCLACHLAN_3_3: return arkodeSymplecticMcLachlan3();
  case ARKODE_SPRK_MCLACHLAN_4_4: return arkodeSymplecticMcLachlan4();
  case ARKODE_SPRK_CANDY_ROZMUS_4_4: return arkodeSymplecticCandyRozmus4();
  case ARKODE_SPRK_MCLACHLAN_5_6: return arkodeSymplecticMcLachlan5();
  case ARKODE_SPRK_YOSHIDA_6_8: return arkodeSymplecticYoshida6();
  case ARKODE_SPRK_SUZUKI_UMENO_8_16: return arkodeSymplecticSuzukiUmeno816();
  case ARKODE_SPRK_SOFRONIOU_10_36: return arkodeSymplecticSofroniou10();
  default: return NULL;
  }
}

ARKodeSPRKTable ARKodeSPRKTable_LoadByName(const char* method)
{
  if (!strcmp(method, "ARKODE_SPRK_EULER_1_1"))
  {
    return arkodeSymplecticEuler();
  }
  if (!strcmp(method, "ARKODE_SPRK_LEAPFROG_2_2"))
  {
    return arkodeSymplecticLeapfrog2();
  }
  if (!strcmp(method, "ARKODE_SPRK_PSEUDO_LEAPFROG_2_2"))
  {
    return arkodeSymplecticPseudoLeapfrog2();
  }
  if (!strcmp(method, "ARKODE_SPRK_RUTH_3_3"))
  {
    return arkodeSymplecticRuth3();
  }
  if (!strcmp(method, "ARKODE_SPRK_MCLACHLAN_2_2"))
  {
    return arkodeSymplecticMcLachlan2();
  }
  if (!strcmp(method, "ARKODE_SPRK_MCLACHLAN_3_3"))
  {
    return arkodeSymplecticMcLachlan3();
  }
  if (!strcmp(method, "ARKODE_SPRK_MCLACHLAN_4_4"))
  {
    return arkodeSymplecticMcLachlan4();
  }
  if (!strcmp(method, "ARKODE_SPRK_CANDY_ROZMUS_4_4"))
  {
    return arkodeSymplecticCandyRozmus4();
  }
  if (!strcmp(method, "ARKODE_SPRK_MCLACHLAN_5_6"))
  {
    return arkodeSymplecticMcLachlan5();
  }
  if (!strcmp(method, "ARKODE_SPRK_YOSHIDA_6_8"))
  {
    return arkodeSymplecticYoshida6();
  }
  if (!strcmp(method, "ARKODE_SPRK_SUZUKI_UMENO_8_16"))
  {
    return arkodeSymplecticSuzukiUmeno816();
  }
  if (!strcmp(method, "ARKODE_SPRK_SOFRONIOU_10_36"))
  {
    return arkodeSymplecticSofroniou10();
  }
  return NULL;
}

ARKodeSPRKTable ARKodeSPRKTable_Copy(ARKodeSPRKTable that_sprk_table)
{
  int i                      = 0;
  ARKodeSPRKTable sprk_table = NULL;

  sprk_table = ARKodeSPRKTable_Alloc(that_sprk_table->stages);

  sprk_table->q = that_sprk_table->q;

  for (i = 0; i < sprk_table->stages; ++i)
  {
    sprk_table->ahat[i] = that_sprk_table->ahat[i];
    sprk_table->a[i]    = that_sprk_table->a[i];
  }

  return sprk_table;
}

void ARKodeSPRKTable_Space(ARKodeSPRKTable sprk_table, sunindextype* liw,
                           sunindextype* lrw)
{
  *liw = 2;
  *lrw = sprk_table->stages * 2;
}

void ARKodeSPRKTable_Free(ARKodeSPRKTable sprk_table)
{
  if (sprk_table)
  {
    if (sprk_table->ahat) { free(sprk_table->ahat); }
    if (sprk_table->a) { free(sprk_table->a); }
    free(sprk_table);
  }
}

void ARKodeSPRKTable_Write(ARKodeSPRKTable sprk_table, FILE* outfile)
{
  ARKodeButcherTable a = NULL;
  ARKodeButcherTable b = NULL;

  ARKodeSPRKTable_ToButcher(sprk_table, &a, &b);

  ARKodeButcherTable_Write(a, outfile);
  ARKodeButcherTable_Write(b, outfile);

  ARKodeButcherTable_Free(a);
  ARKodeButcherTable_Free(b);
}

int ARKodeSPRKTable_ToButcher(ARKodeSPRKTable sprk_table,
                              ARKodeButcherTable* a_ptr,
                              ARKodeButcherTable* b_ptr)
{
  int i                = 0;
  int j                = 0;
  ARKodeButcherTable a = NULL;
  ARKodeButcherTable b = NULL;

  a = ARKodeButcherTable_Alloc(sprk_table->stages, SUNFALSE);
  if (!a) { return ARK_MEM_FAIL; }
  b = ARKodeButcherTable_Alloc(sprk_table->stages, SUNFALSE);
  if (!b)
  {
    if (a) { ARKodeButcherTable_Free(a); }
    return ARK_MEM_FAIL;
  }

  /* DIRK table */
  for (i = 0; i < sprk_table->stages; ++i)
  {
    b->b[i] = sprk_table->ahat[i];
    for (j = 0; j <= i; ++j) { b->A[i][j] = sprk_table->ahat[j]; }
    /* Time weights: C_j = sum_{i=0}^{j} b_i */

    /* Time weights: C_j = sum_{i=0}^{j-1} b_i */
    for (j = 0; j < sprk_table->stages; ++j)
    {
      for (i = 0; i <= j; ++i) { b->c[j] += sprk_table->ahat[i]; }
    }

    /* Explicit table */
    for (i = 0; i < sprk_table->stages; ++i)
    {
      a->b[i] = sprk_table->a[i];
      for (j = 0; j < i; ++j) { a->A[i][j] = sprk_table->a[j]; }
    }

    /* Time weights: c_j = sum_{i=0}^{j-1} a_i */
    for (j = 0; j < sprk_table->stages; ++j)
    {
      for (i = 0; i < j; ++i) { a->c[j] += sprk_table->a[i]; }
    }

    /* Set method order */
    a->q = sprk_table->q;
    b->q = sprk_table->q;

    /* No embedding, so set embedding order to 0 */
    a->p = 0;
    b->p = 0;
  }

  *a_ptr = a;
  *b_ptr = b;

  return ARK_SUCCESS;
}
