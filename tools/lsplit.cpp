#include <cctype>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

std::string str_tolower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); }
                   );
    return s;
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " (LAPACK file with dependencies)\n";
        return 1;
    }

    std::fstream in;
    in.open(argv[1]);

    if (!in.is_open()) {
        std::cerr << "Could not open file " << argv[1] << " for reading\n";
        return 2;
    }

    static constexpr int BUFSIZE = 1024;
    char buf[BUFSIZE];
    std::string text;
    std::fstream out;
    
    while (in.getline(buf, BUFSIZE)) {
        text = buf;
        
        // look for start of next file marker
        auto n = text.find("*> \\brief \\b ");
        if (n == std::string::npos) {
            if (out.is_open()) out << text << "\n";
        } else {
            if (out.is_open()) out.close();
            auto e = text.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890", 13);
            std::string name;
            if (e == std::string::npos)
                name = text.substr(13, e);
            else
                name = text.substr(13, e - 13);
            name = str_tolower(name);
            name += ".f";
            out.open(name, std::ios_base::out);
            std::cout << "Found: " << name << "\n";
            out << text << "\n";
        }
    }
    std::cout << "Done\n";
    if (out.is_open()) out.close();
    in.close();
}        

    
