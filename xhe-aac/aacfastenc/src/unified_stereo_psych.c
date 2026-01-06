
/* -----------------------------------------------------------------------------
Software Copyright License for The Fraunhofer FDK Extended High Efficiency AAC
Encoder Software for Android

© Copyright 1995 - 2025 Fraunhofer-Gesellschaft zur Förderung der angewandten
Forschung e.V. and Contributors
All rights reserved.

1.    INTRODUCTION

The Fraunhofer FDK Extended High Efficiency AAC Encoder Software for Android
("FDK Extended High Efficiency AAC Encoder") is software that implements the
encoding of digital audio according to the MPEG-D Unified Speech and Audio
Coding (USAC) standard and MPEG-D Dynamic Range Control (DRC) standard. This FDK
Extended High Efficiency AAC Encoder Software is intended to be used on a wide
variety of Android devices. It is technically not suited to encode content for
digital radio broadcasting services, including DRM and similar standards.

Patent licenses for necessary patent claims for the FDK Extended High Efficiency
AAC Encoder Software (including those of Fraunhofer), for the use in commercial
products and services, may be obtained from the respective patent owners
individually and/or from Via Licensing Alliance (www.via-la.com).

Fraunhofer supports the development of Extended High Efficiency AAC products and
services by offering additional software, documentation, and technical advice.
In addition, it operates the xHE-AAC Trademark Program to ease interoperability
testing of end products. Please visit http://www.xhe-aac.com for more
information.

2.    COPYRIGHT LICENSE

Redistribution and use in source and binary forms, with or without modification,
are permitted without payment of copyright license fees, provided that you
satisfy the following conditions:

You must retain the complete text of this software license in redistributions of
the FDK Extended High Efficiency AAC Encoder Software or your modifications
thereto in source code form.

You must retain the complete text of this software license in the documentation
and/or other materials provided with redistributions of the FDK Extended High
Efficiency AAC Encoder Software or your modifications thereto in binary form.
You must make available free of charge copies of the complete source code of the
FDK Extended High Efficiency AAC Encoder Software and your modifications thereto
to recipients of copies in binary form.

The name of Fraunhofer may not be used to endorse or promote products derived
from this software without prior written permission.

You may not charge copyright license fees for anyone to use, copy or distribute
the FDK Extended High Efficiency AAC Encoder Software or your modifications
thereto.

Your modified versions of the FDK Extended High Efficiency AAC Encoder Software
must carry prominent notices stating that you changed the software and the date
of any change. For modified versions of the FDK Extended High Efficiency AAC
Encoder Software, the term "Fraunhofer FDK Extended High Efficiency AAC Encoder
Software for Android" must be replaced by the term "Third-Party Modified Version
of the Fraunhofer FDK Extended High Efficiency AAC Encoder Software for
Android."

3.    NO PATENT LICENSE

NO EXPRESS OR IMPLIED LICENSES TO ANY PATENT CLAIMS, including without
limitation the patents of Fraunhofer, ARE GRANTED BY THIS SOFTWARE LICENSE.
Fraunhofer provides no warranty for patent non-infringement with respect to this
software. You may use this FDK Extended High Efficiency AAC Encoder Software or
modifications thereto only for purposes that are authorized by appropriate
patent licenses.

4.    DISCLAIMER

This FDK Extended High Efficiency AAC Encoder Software is provided by Fraunhofer
on behalf of the copyright holders and contributors "AS IS" and WITHOUT ANY
EXPRESS OR IMPLIED WARRANTIES, including but not limited to the implied
warranties of merchantability and fitness for a particular purpose. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE for any direct, indirect,
incidental, special, exemplary, or consequential damages, including but not
limited to procurement of substitute goods or services; loss of use, data, or
profits, or business interruption, however caused and on any theory of
liability, whether in contract, strict liability, or tort (including
negligence), arising in any way out of the use of this software, even if advised
of the possibility of such damage.

5.    CONTACT INFORMATION

Fraunhofer Institute for Integrated Circuits IIS
Attention: Division Audio and Media Technologies - FDK Extended High Efficiency
AAC Encoder
Am Wolfsmantel 33
91058 Erlangen, Germany

www.iis.fraunhofer.de/amm
amm-info@iis.fraunhofer.de
----------------------------------------------------------------------------- */

#include <float.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "mathlib.h"
#include "glob_con.h"
#include "psy_const.h"
#include "unified_stereo_psych.h"
#include "interface.h"

void iisaacfenc_UniStePsyProcessing(
    float *sfbThresholdLeft,
    float *sfbThresholdRight,
    int *msDigest,
    int *msMask,
    const int sfbCnt,
    const int sfbPerGroup,
    const int maxSfbPerGroup,
    const int dualMono,
    const float *invPredGain) {
  int sfb, sfboffs;

  *msDigest = MS_NONE;
  setINT(0, msMask, sfbCnt);

  if (dualMono) {
    return;
  }

  for (sfboffs = 0; sfboffs < sfbCnt; sfboffs += sfbPerGroup) {
    for (sfb = sfboffs; sfb < sfboffs + maxSfbPerGroup; sfb++) {
      float energyCorrFac = 1.0f / (1.0f + invPredGain[sfb]);
      float thrUni = energyCorrFac * (invPredGain[sfb] * sfbThresholdLeft[sfb] + sfbThresholdRight[sfb]);
      sfbThresholdLeft[sfb] = max(sfbThresholdLeft[sfb], thrUni);
      sfbThresholdRight[sfb] = max(sfbThresholdRight[sfb], thrUni);
    }
  }
}

void iisaacfenc_mapUniSteCldToPredGain(
    const HANDLE_UNISTE hUniSte,
    const PSY_CONFIGURATION *psyConf,
    const int windowSequence,
    const int sfbCnt,
    const int sfbPerGroup,
    const int maxSfbPerGroup,
    const int *groupLen,
    float uniStePredGainSfb[MAX_SFB]) {
  int i;
  int k;
  float paramBandBorders[28 + 1];
  int sfb;
  int sfboffs;
  float umxMatSfbRe[MAX_SFB][4];
  float umxMatSfbIm[MAX_SFB][4];

  if (hUniSte == 0) {
    for (sfb = 0; sfb < sfbCnt; sfb += sfbPerGroup) {
      for (sfboffs = 0; sfboffs < maxSfbPerGroup; sfboffs++) {
        uniStePredGainSfb[sfb + sfboffs] = 1.0f;
      }
    }
  } else {
    int timeSlot;
    static const int LS_TRANS = ((FRAME_LEN_LONG - FRAME_LEN_SHORT) / 2);
    int pos;
    int groupAcc = 0;

    if (windowSequence == SHORT_WINDOW) {
      memcpy(paramBandBorders, hUniSte->paramBandBordersShort, (hUniSte->nParamBands + 1) * sizeof(float));
    } else {
      memcpy(paramBandBorders, hUniSte->paramBandBordersLong, (hUniSte->nParamBands + 1) * sizeof(float));
    }

    switch (windowSequence) {
      case LONG_WINDOW:
      case STOPSTART_WINDOW:
        timeSlot = (2 * FRAME_LEN_LONG) / 64;
        break;
      case START_WINDOW:
      case STOP_WINDOW:
        timeSlot = (FRAME_LEN_LONG + LS_TRANS) / 64;
        break;
      default:
        timeSlot = -1;
        break;
    }

    for (sfb = 0; sfb < sfbCnt; sfb += sfbPerGroup) {
      if (windowSequence == SHORT_WINDOW) {
        timeSlot = (FRAME_LEN_LONG + 64 + groupAcc * FRAME_LEN_SHORT) / 64;
      }

      for (sfboffs = 0; sfboffs < maxSfbPerGroup; sfboffs++) {
        int sfbStartLine = psyConf->sfbOffset[sfboffs];
        int sfbEndLine = min(psyConf->lowpassLine, psyConf->sfbOffset[sfboffs + 1]);
        int paramBandFirst = -1;
        int paramBandLast = -1;
        assert(sfb + sfboffs < MAX_SFB);

        for (i = 0; i < hUniSte->nParamBands; i++) {
          if (paramBandBorders[i] > sfbStartLine) {
            paramBandFirst = i - 1;
            break;
          }
        }

        for (i = 0; i < hUniSte->nParamBands; i++) {
          if (paramBandBorders[i] >= sfbEndLine) {
            paramBandLast = i - 1;
            break;
          }
        }

        if ((paramBandFirst < hUniSte->nResidualBands) &&
            (paramBandLast >= hUniSte->nResidualBands)) {
          paramBandBorders[paramBandLast] = (float)sfbEndLine;
          paramBandLast -= 1;
        }

        uniStePredGainSfb[sfb + sfboffs] = 0.0f;
        for (k = 0; k < 4; k++) {
          umxMatSfbRe[sfb + sfboffs][k] = umxMatSfbIm[sfb + sfboffs][k] = 0.0f;
        }

        pos = sfbStartLine;
        for (i = paramBandFirst; i <= paramBandLast; i++) {
          for (k = 0; k < 4; k++) {
            umxMatSfbRe[sfb + sfboffs][k] += (min(sfbEndLine, paramBandBorders[i + 1]) - pos) * hUniSte->umxMatRe[timeSlot][i][k] / (sfbEndLine - sfbStartLine);
            umxMatSfbIm[sfb + sfboffs][k] += (min(sfbEndLine, paramBandBorders[i + 1]) - pos) * hUniSte->umxMatIm[timeSlot][i][k] / (sfbEndLine - sfbStartLine);
          }

          uniStePredGainSfb[sfb + sfboffs] += (float)pow(10, (hUniSte->cld[timeSlot][i] > 0 ? -hUniSte->cld[timeSlot][i] : hUniSte->cld[timeSlot][i]) * 0.1f) * (min(sfbEndLine, paramBandBorders[i + 1]) - pos) / (sfbEndLine - sfbStartLine);

          pos = (int)paramBandBorders[i + 1];
        }

        uniStePredGainSfb[sfb + sfboffs] = min(
            (umxMatSfbRe[sfb + sfboffs][0] * umxMatSfbRe[sfb + sfboffs][0] +
             umxMatSfbIm[sfb + sfboffs][0] * umxMatSfbIm[sfb + sfboffs][0]),
            (umxMatSfbRe[sfb + sfboffs][2] * umxMatSfbRe[sfb + sfboffs][2] +
             umxMatSfbIm[sfb + sfboffs][2] * umxMatSfbIm[sfb + sfboffs][2]));
      }
      groupAcc += groupLen[sfb / sfbPerGroup];
    }
  }
}
