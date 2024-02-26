//
// Created by Roman Wallner- Silberhuber on 21.02.24.
//
#import "SuperElement.h"

SuperElement::SuperElement(int size) : arraySize(size)
{
    members = new int[arraySize];
}

SuperElement::SuperElement()
= default;

int SuperElement::getMember(int index) const
{
    if (usingDynamic)
    {
        return membersVector[index];
    }
    else
    {
        return members[index];
    }
}

void SuperElement::setMember(int index, int value)
{
    if (usingDynamic)
    {
        membersVector[index] = value;
    }
    else
    {
        members[index] = value;
    }
}

void SuperElement::addMember(int value)
{
    if (!usingDynamic)
    {
        // Handle potential error if you intend to limit to fixed-size
        return;
    }
    membersVector.push_back(value);
}

SuperElement::~SuperElement()
{
    if (!usingDynamic)
    {
        delete[] members; // Deallocate fixed-size array
    }
}