#ifndef GRAPHICSCONTROLLER_H
#define GRAPHICSCONTROLLER_H

#include "NURBSVisualization.h"
#include "raylib.h"
#include <tinynurbs/tinynurbs.h>

namespace GraphicsController
{

void
InitializeWindow (int width, int height, const char *title);
void
RenderLoop (const tinynurbs::Curve<float> &crv);
void
Terminate ();

}

#endif // GRAPHICSCONTROLLER_H
