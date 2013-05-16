#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <fstream>

class Vector2D {
  public:
    Vector2D(): x_(0), y_(0) {}
    void x(double x) {x_ = x;}
    void y(double y) {y_ = y;}
    double x() { return x_; }
    double y() { return y_; }
    Vector2D operator-(Vector2D rhs) {
        Vector2D ret;
        ret.x( x_ - rhs.x() );
        ret.y( y_ - rhs.y() );
        return ret;
    }
    Vector2D operator+=(Vector2D& rhs) {
        x_ += rhs.x();
        y_ += rhs.y();
        return *this;
    }
    Vector2D operator*=(double rhs) {
        x_ *= rhs;
        y_ *= rhs;
        return *this;
    }
    double distance() {
        return sqrt(x_ * x_ + y_ * y_);
    }
    void normalize() {
        double magnitude = sqrt(x_ * x_ + y_ + y_);
        x_ /= magnitude;
        y_ /= magnitude;
    }
  private:
    double x_, y_;
};

struct part {
    double x, y, vx, vy, carga, raio;
};

typedef struct part particula;

int main(int argc, char* argv[]) {
    std::ifstream file(argv[1], std::ifstream::in);
    if(!file) return 1;
    int n_partic = 0;
    double size_x = 0, size_y = 0;
    double c = atof(argv[2]);
    double epsilon = atof(argv[3]);
    double tau = atof(argv[4]);
    int iters = atoi(argv[5]);
    int procs = atoi(argv[6]);
    int i, j, k;
    
    file >> n_partic;
    file >> size_x;
    file >> size_y;
    particula particulas[n_partic];
    Vector2D forcas[n_partic];
    for(i = 0; i < n_partic; ++i) {
        file >> particulas[i].x; file >> particulas[i].y;
        file >> particulas[i].vx; file >> particulas[i].vy;
        file >> particulas[i].carga; file >> particulas[i].raio;
    }
    printf("%s %f %f %f %d %d\n", argv[1], c, epsilon, tau, iters, procs);
    printf("%d, %f, %f\n", n_partic, size_x, size_y);
    for(i = 0; i < n_partic; ++i)
        printf("\t%f %f %f %f %f %f\n", particulas[i].x, particulas[i].y, particulas[i].vx, particulas[i].vy, particulas[i].carga,particulas[i].raio);
    
    #pragma omp parallel
    for(i = 0; i < iters; ++i) {
        #pragma omp for
        for(j = 0; j < n_partic; ++j) {
            forcas[j].x(0);
            forcas[j].y(0);
            printf("%d\n", j);
            for(k = 0; k < n_partic; ++k)
                if(j == k)
                    continue;
                Vector2D vec;
                vec.x(particulas[j].x - particulas[k].y);
                vec.y(particulas[j].y - particulas[k].y);
                double distance = vec.distance();
                if(distance > epsilon)
                    continue;
                forcas[j] += (vec * distance);
        }
    }
}