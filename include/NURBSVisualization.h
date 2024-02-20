//
// Created by Roman Wallner- Silberhuber on 17.02.24.
//

#ifndef _NURBSVISUALIZATION_H_
#define _NURBSVISUALIZATION_H_

#include "tinynurbs/tinynurbs.h"
#include "raylib.h"

namespace NURBSVisualization {

void DrawCurve(const tinynurbs::Curve<float>& crv);

}
#endif //_NURBSVISUALIZATION_H_
