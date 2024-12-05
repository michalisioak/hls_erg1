#ifndef __READ_TXT_H__
#define __READ_TXT_H__
#include "ac_int.h"
#include <fstream>
#include <iostream>
#include <sstream>

template <typename Type, int N, int L>
void read_data(Type mat[N][L], std::string filename) {
  std::ifstream myFile(filename);
  if (myFile.is_open()) { // throw std::runtime_error("Could not open file"); //
                          // Make sure the file is open

    std::string line;
    int val;
    int rowIdx = 0;

    // Read data, line by line
    while (std::getline(myFile, line)) {
      std::stringstream ss(line); // Create a stringstream of the current line
      int colIdx = 0;             // Keep track of the current column index

      // Extract each integer
      while (ss >> val) {
        mat[rowIdx][colIdx] = val; // Write current input value
        if (ss.peek() == ',')
          ss.ignore(); // If the next token is a comma, ignore it and move on

        colIdx++; // Increment the Column index
        if (colIdx == L)
          break;
      }
      rowIdx++; // Increment the Row index
      if (rowIdx == N)
        break;
    }
    myFile.close(); // Close file
  }
}

#endif