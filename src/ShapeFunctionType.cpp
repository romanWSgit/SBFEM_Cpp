//
// Created by Roman Wallner- Silberhuber on 20.02.24.
//
#include "ShapeFunctionType.h"



std::string shapeFunctionTypeToString(ShapeFunctionType type)
{
    switch (type)
    {
    case ShapeFunctionType::STANDARD:
        return "standard shape functions";
    case ShapeFunctionType::HIERARCHICAL:
        return "hierarchical shape functions";
    default:
        return "unknown shape function type";
    }
}