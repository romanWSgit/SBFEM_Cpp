#include "GraphicsController.h"
// Include your NURBS visualization header

namespace GraphicsController
{
void
InitializeWindow (int width, int height, const char *title)
{
    InitWindow (width, height, title);
    SetTargetFPS (60);
}

void
RenderLoop (const tinynurbs::Curve<float> &crv)
{
    while (!WindowShouldClose ())
        {
            BeginDrawing ();
            ClearBackground (RAYWHITE);

            // Draw your curve here using the provided NURBS curve
            NURBSVisualization::DrawCurve (crv);

            EndDrawing ();
        }
}

void
Terminate ()
{
    CloseWindow ();
}

}