#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

vector<vector<double>> loadMTXFile(string filename){
    
    ifstream infile(filename);
    vector<vector<double>> mat;

    if(infile.is_open()){
        string line;
        bool matSize = false;

        while(getline(infile, line)){
            
            if(line.at(0) == '%')
                continue;

            istringstream iss(line);
            if(!matSize){
                int rowNum, colNum, total;
                iss >> rowNum >> colNum >> total;
                mat.reserve(rowNum);

                for(int i = 0; i < rowNum; i++){
                    mat.push_back(vector<double>(colNum));
                }
                matSize = true;
            }
            else{
                int r, c;
                double entry;
                iss >> r >> c >> entry;
                mat[r-1][c-1] = entry;
            }
        }
    }
    else {
        cerr << "Failed to open file: " << filename << endl;
    }
    return mat;
}

void writeMTXFile(string filename, vector<vector<double>> mat){
    
    ofstream outfile(filename);

    int rowNum = mat.size();
    int colNum = mat[0].size();

    if(outfile.is_open()){
        outfile << rowNum << " " << colNum << " " << rowNum*colNum << endl;

        for(int i = 1; i <= rowNum; i++){
            for(int j = 1; j <= colNum; j++){ 
                outfile << i << " " << j << " " << mat[i-1][j-1] << endl;
            }
        }
    }
    else{
        cerr << "Failed to open file: " << filename << endl;
    }
}




