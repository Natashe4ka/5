#include "head.hpp"

int main() {
    Parametrs param;
    std::vector<Figure*> objects; 

    std::string file;

    std::cout << "Enter file with parametrs: ";
    std::cin >> file;
    
    bool res1 = GetParametrs(file, param);
    if(!res1) {
        return -1;
    }

    std::cout << "Enter file with objects: ";
    std::cin >> file;
    
    bool res2 = GetFigures(file, objects);
    if(!res2) {
        return -1;
    }

    
    std::map<double, size_t> cntVolumeFig;

    for(size_t i = 0; i < objects.size(); ++i) {
        double v = objects[i]->calculate_volume();
        ++cntVolumeFig[v];
    }
    
    Vec3d colour_min(64, 64, 64);
    Vec3d colour_max(191, 191, 191);

    double b = (191-64)/(cntVolumeFig.size()-1);
    for(size_t i = 0; i < objects.size(); ++i) {
        double q = 0; // порядковый номер
        for(auto it : cntVolumeFig) {
            if(it.first == objects[i]->calculate_volume()) break;
            ++q;
        }
        objects[i]->UpdateColour(colour_min+(Vec3d(b, b, b)*q));
    }

    create_image(objects, param);

    return 0;
}

