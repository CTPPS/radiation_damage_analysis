//CSVParser.h
//****************************************************
//This lib contains the exception class for personal 
//error handling
//****************************************************

#ifndef MYEXCEPTION_H
#define MYEXCEPTION_H

#include <string>

struct MyException : public std::exception
{
   std::string s;

   MyException(std::string ss) : s(ss) {}
   ~MyException() throw () {}

   const char* what() const throw()
   {
	   return s.c_str();
   }
};

#endif
