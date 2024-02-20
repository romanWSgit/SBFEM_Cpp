//
// Created by Roman Wallner- Silberhuber on 20.02.24.
//

#ifndef SHAPEFUNCTIONTYPE_H_
#define SHAPEFUNCTIONTYPE_H_

#include <iostream>
#include <string>

enum class ShapeFunctionType
{
    STANDARD,
    HIERARCHICAL
};

std::string shapeFunctionTypeToString(ShapeFunctionType type);


#endif //SHAPEFUNCTIONTYPE_H_
