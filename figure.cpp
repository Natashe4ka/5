#include "Figure.hpp"


Vec3d Figure::colour() const {
    return _colour;
};
void Figure::UpdateColour(const Vec3d& colour) {
    _colour = colour*(1/255.0);
}
double Figure::volume() const {
    return _volume;
};


double Sphere::ray_intersect(const Vec3d& orig, const Vec3d& dir) const {
    Vec3d L = orig - _center;

    double a = 1.0; // (dir; dir) = 1
    double b = 2*(dir*L);
    double c = (L*L) - _R*_R;
        
    double det = b*b-4*a*c;

    if(det > 0) {
        double t1 = (-b + sqrt(det)) / 2;
        double t2 = (-b - sqrt(det)) / 2;
        if(t1 > 0 && t2 > 0) {
            return std::min(t1, t2);
        }
    }
        
    return -1;
}
Vec3d Sphere::center() const {
    return _center;
}
double Sphere::R() const {
    return _R;
}
double Sphere::calculate_volume() {
    _volume = (4/3.0)*PI*_R*_R*_R;
    return _volume;
}



double Box::ray_intersect(const Vec3d& orig, const Vec3d& dir) const {
    double t_near = std::numeric_limits<double>::min(),
            t_far  = std::numeric_limits<double>::max();
    double t1, t2;

    for(int i = 0; i < 3; ++i) {
        double d = dir[i];
        if(d < 0.0) d *= -1;
        if(d >= std::numeric_limits<double>::epsilon()) {
            t1 = (_A[i] - orig[i]) / dir[i];   
            t2 = (_B[i] - orig[i]) / dir[i];    
            if(t1 > t2) std::swap(t1, t2);


            if(t1 > t_near) t_near = t1;
            if(t2 < t_far) t_far = t2;

            if(t_near > t_far) return -1;
            if(t_far < 0.0) return -1;
        }
        else {
            if(orig[i] < _A[i] || orig[i] > _B[i] ) return -1;
        }
    }
    //std::cout << t_near << "\n";
    if(t_near <= t_far && t_far >=0) return t_near;
    else return -1;
}
Vec3d Box::A() const {
    return _A;
}
Vec3d Box::B() const {
    return _B;
}
Vec3d Box::center() const {
    return Vec3d( (_A[0]+_B[0])/2, (_A[1]+_B[1])/2, (_A[2]+_B[2])/2);
}
double Box::calculate_volume() {
    double a = _B[0]-_A[0];
    double b = _B[1]-_A[1];
    double c = _B[2]-_A[2];

    if(a < 0.0) a = -a;
    if(b < 0.0) b = -b;
    if(c < 0.0) c = -c;

    _volume = a*b*c;
    return _volume;
}


double Tetrahedron::ray_intersect(const Vec3d& orig, const Vec3d& dir) const {
        double t1 = triangle_intersection(orig, dir, _A, _B, _C);
        double t2 = triangle_intersection(orig, dir, _A, _B, _D);
        double t3 = triangle_intersection(orig, dir, _D, _B, _C);
        double t4 = triangle_intersection(orig, dir, _A, _D, _C);

        double t_leave = std::max(t1, std::max(t2, std::max(t3, t4)));

        if(t_leave <= 0.0) return -1;

        double t_enter = -1.0;

        if(t1 != 0.0) t_enter = t1;
        if(t2 != 0.0) t_enter = std::min(t_enter, t2);
        if(t3 != 0.0) t_enter = std::min(t_enter, t3);
        if(t4 != 0.0) t_enter = std::min(t_enter, t4);

        return t_enter;
}
Vec3d Tetrahedron::center() const {
    double sumi[3] = {0.0, 0.0, 0.0};
    for(size_t i = 0; i < 3; ++i) {
        sumi[i] += (_A[i]+_B[i]+_C[i]+_D[i]);
    }
    return Vec3d(sumi[0]/4, sumi[1]/4, sumi[2]/4);
}
double Tetrahedron::calculate_volume() {
    double V = ((_A-_D)*cross(_B-_D, _C-_D))*1/6;
    _volume = V >= 0.0 ? V : -V;
    return _volume;
}
