#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include "graph.hpp"
int main(int argc, char *argv[])
{
    const char *infile = argv[1];
    const char *outfile = argv[2];
    if (argc < 3)
    {
        std::cout << "format: ./treat_datafile \"inputfile\" \"outputfile\"\n";
        return -1;
    }
    ifstream weightFile("input/weightFile1.txt");
    std::ifstream input(infile);
    assert(input.is_open());
    //open output data file

    std::ofstream output(outfile);
    output.clear();

    std::string line;

    // assert(input.is_open());
    srand(time(0));
    int N;
    for (int i = 0; i < 17; i++)
    {
        //read the header
        getline(input, line);
        if (i == 4)
        {
            std::stringstream s(line);
            std::string word;
            s >> word;
            s >> word;
            s >> word;
            std::stringstream(word) >> N;
        }
    }
    std::cout << N << "\n";
    output << N << '\n';
    graph g(N, &input, &weightFile);
       
        std::cout <<"oi\n";

    g.graphToFile(&output);

    output.close();
    return 0;
}