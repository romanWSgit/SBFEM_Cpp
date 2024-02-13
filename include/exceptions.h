//
// Created by Roman Wallner- Silberhuber on 18.07.23.
//

#ifndef SBFEM_EXCEPTIONS_H
#define SBFEM_EXCEPTIONS_H

#include <exception>
#include <string>

class MyBaseException : public std::exception
{
    // Common code for all exceptions goes here
    std::string message;
public:
    MyBaseException(const std::string& msg) : message(msg)
    {
    }
    const char* what() const noexcept override {
        return message.c_str();
    }
};

class MyFileNotFoundException : public MyBaseException
{
public:
    MyFileNotFoundException(const std::string& msg) : MyBaseException(msg)
    {
    }
};

class MyNetworkException : public MyBaseException
{
public:
    MyNetworkException(const std::string& msg) : MyBaseException(msg) {}
};


#endif //SBFEM_EXCEPTIONS_H
