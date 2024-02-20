//
// Created by Roman Wallner- Silberhuber on 17.02.24.
//
#include "NURBSVisualization.h"

void
NURBSVisualization::DrawCurve (const tinynurbs::Curve<float> &crv)
{
    float scale = 200.0f; // Adjust as necessary for your visualization
    Vector2 windowCenter
        = { GetScreenWidth () / 2.0f, GetScreenHeight () / 2.0f };
    float lineThickness = 3.0f; // Increase this value for a bolder curve

    // Draw the curve
    for (int i = 0; i < 100; ++i)
        {
            float t1 = crv.knots[crv.degree]
                       + (i / 100.0f)
                             * (crv.knots[crv.knots.size () - crv.degree - 1]
                                - crv.knots[crv.degree]);
            float t2 = crv.knots[crv.degree]
                       + ((i + 1) / 100.0f)
                             * (crv.knots[crv.knots.size () - crv.degree - 1]
                                - crv.knots[crv.degree]);

            glm::vec3 pt1_3D
                = tinynurbs::curvePoint (crv, t1); // Evaluate curve at t1
            glm::vec3 pt2_3D
                = tinynurbs::curvePoint (crv, t2); // Evaluate curve at t2

//            Vector2 start = { windowCenter.x + pt1_3D.x * scale,
//                              windowCenter.y + pt1_3D.y * scale };
//            Vector2 end = { windowCenter.x + pt2_3D.x * scale,
//                            windowCenter.y + pt2_3D.y * scale };

            // Adjust Y coordinate by inverting the value relative to the window center
            Vector2 start = {windowCenter.x + pt1_3D.x * scale, windowCenter.y - pt1_3D.y * scale}; // Note the minus sign
            Vector2 end = {windowCenter.x + pt2_3D.x * scale, windowCenter.y - pt2_3D.y * scale}; // Note the minus sign


            //DrawLineV (start, end, BLACK); // Draw curve segment
            DrawLineEx(start, end, lineThickness, BLACK); // Draw curve segment with specified thickness
            // Optionally, draw sample point for t1 (not for t2 to avoid
            // duplicate draws)
//            DrawCircle (start.x, start.y, 3,
//                        RED); // Sample point as a small circle
        }

    // Draw control points
    for (const auto &cp : crv.control_points)
        {
            Vector2 cpPos = { windowCenter.x + cp.x * scale,
                              windowCenter.y + cp.y * scale };
            DrawCircle (cpPos.x, cpPos.y, 5,
                        GREEN); // Control point as a slightly larger circle
            DrawCircleLines (
                cpPos.x, cpPos.y, 5,
                DARKGREEN); // Optionally, add an outline for clarity
        }
}
