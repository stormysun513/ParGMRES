#include <iostream>

void run(char *filename);

int main(int argc, char *argv[]) {
    
    if(argc != 2){
        std::cout << "\nUsage:\n";
        std::cout << "      ./gmres <filename>\n";
        return -1;
    }
    
    run(argv[1]);
    return 0;
}
