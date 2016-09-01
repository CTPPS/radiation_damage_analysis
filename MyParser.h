//CSVParser.h
//****************************************************
//This lib contains the method used to parse multiple 
//lines files and transfer all the lines in a vector 
//of strings
//****************************************************

#ifndef MYPARSER_H
#define MYPARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "MyException.h"

using namespace std;

// TODO: remove
void LineParser(ifstream&, vector<string>&);
void CSVLinesParser(vector<string>, vector<vector<string> >&);

/// Function to turn a multiple lines file into a vector of strings 
void LineParser(ifstream& file_to_parse, vector<string>& lines)
{
  try
  { 
    if (!file_to_parse)
	{
      throw MyException("ERROR: file to parse not found!");
      goto END;     
    }
    
    while (true)
	{
      if (file_to_parse.eof())
	  	break;
      string line;
      getline(file_to_parse,line);
      lines.push_back(line);
    }
    
    END:;
  }

  catch(MyException& caught)
  {
    std::cout << caught.what() << endl;
  }
  
}

/// Function to turn a vector of lines coming from a .csv in a vector
/// of vectors containing the lables
void CSVLinesParser(vector<string> lines, vector<vector<string> >& spreadsheet)
{
  for (vector<string>::iterator line=lines.begin();line!=lines.end();++line)
  {
    string field_content;
    istringstream linestream(*line);
    if (!linestream.eof()){
      spreadsheet.push_back(vector<string>());
      while(1){
	getline(linestream,field_content,',');
	spreadsheet.back().push_back(field_content);
	if (linestream.eof()){break;}
      }
    }
  }
}

#endif
