#ifndef IMAGE_WRITER_H
#define IMAGE_WRITER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

class ImageWriter {
public:
    static void writePPM(const std::string& filename, int width, int height, const std::vector<unsigned char>& data) {
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }

        // PPM Header
        ofs << "P6\n" << width << " " << height << "\n255\n";

        // Data
        ofs.write(reinterpret_cast<const char*>(data.data()), data.size());
        
        ofs.close();
        std::cout << "Image saved to " << filename << std::endl;
    }
};

#endif
