//
// Created by Roman Wallner- Silberhuber on 21.02.24.
//

#ifndef SUPERELEMENT_H_
#define SUPERELEMENT_H_

#include <vector> // For dynamic array option

class SuperElement
{
  public:
    // Constructor for fixed-size array
    SuperElement(int size);

    SuperElement();

    // Accessor (for both fixed and dynamic)
    int getMember(int index) const;

    // Mutator (for both fixed and dynamic)
    void setMember(int index, int value);

    // Add a member to the dynamic array
    void addMember(int value);

    // Destructor (IMPORTANT!)
    ~SuperElement();

  private:
    int arraySize{};
    int *members{};                   // For fixed-size array
    std::vector<int> membersVector; // For dynamic array
    bool usingDynamic = false;      // Flag to track which storage is used
};

#endif // SUPERELEMENT_H_
