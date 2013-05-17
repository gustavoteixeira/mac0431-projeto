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
    void add(double x, double y) {
        x_ += x;
        y_ += y;
    }
  private:
    double x_, y_;
};

struct part {
    double x, y, vx, vy, carga, raio;
};

typedef struct part particula;

int main(int argc, char* argv[]) {
    double x, y, distancia_sqr, distancia, forca = 0;
    int i, j, k, tid;
    std::ifstream file(argv[1], std::ifstream::in);
    if(!file) return 1;
    int n_partic = 0;
    double size_x = 0, size_y = 0;
    double c = atof(argv[2]);
    double epsilon = atof(argv[3]);
    double tau = atof(argv[4]);
    int iters = atoi(argv[5]);
    int procs = atoi(argv[6]);
    
    if (argc > 7)
      omp_set_num_threads(atoi(argv[6]));

    file >> n_partic;
    file >> size_x;
    file >> size_y;
    particula particulas[n_partic];
    Vector2D forcas[n_partic];
    int iterador[n_partic];
    double lx[n_partic], ly[n_partic];
    for(i = 0; i < n_partic; ++i) {
        file >> particulas[i].x; file >> particulas[i].y;
        file >> lx[i]; file >> ly[i];
        file >> particulas[i].carga; file >> particulas[i].raio;
    }
    printf("%s %f %f %f %d %d\n", argv[1], c, epsilon, tau, iters, procs);
    printf("%d, %f, %f\n", n_partic, size_x, size_y);
    for(i = 0; i < n_partic; ++i) {
        printf("\t%f %f %f %f %f %f\n", particulas[i].x, particulas[i].y, particulas[i].vx, particulas[i].vy, particulas[i].carga,particulas[i].raio);
        lx[i] = 0; ly[i] = 0;
    }
    
    #pragma omp parallel private(lx, ly, tid)
    {
        tid=omp_get_thread_num();
        for(i = 0; i < iters; ++i) {
            if(i%25000 == 0 || i > 2499900)
                printf("Iteracao %d\n", i);
            #pragma omp for private(x, y, distancia, distancia_sqr, k)
            for(j = 0; j < n_partic; ++j) {
                forcas[j].x(0);
                forcas[j].y(0);
                particulas[j].vx = lx[j];
                particulas[j].vy = ly[j];
                for(k = 0; k < n_partic; ++k) {
                    //#pragma omp critical
                    //printf("-%d, %d-\n", j, k);
                    if(j == k)
                        continue;
                    x = (particulas[j].x - particulas[k].x);
                    y = (particulas[j].y - particulas[k].y);
                    distancia_sqr = x*x + y*y;
                    distancia = sqrt(distancia_sqr);
                    if(distancia > epsilon)
                        continue;
                    forca = c/distancia_sqr;
                    forcas[j].add((x/distancia_sqr)*forca, (y/distancia_sqr)*forca);
                    lx[j] = particulas[j].vx + forcas[j].x();
                    ly[j] = particulas[j].vy + forcas[j].y();
                }
            }
            
            #pragma omp barrier
            #pragma omp for
            for(j = 0; j < n_partic; ++j) {
                particulas[j].x += particulas[j].vx;
                particulas[j].y += particulas[j].vy;
                if(particulas[j].x > size_x)
                    particulas[j].x = 2*size_x - particulas[j].x;
                else if(particulas[j].x < 0)
                    particulas[j].x *= -1;
                
                if(particulas[j].y > size_y)
                    particulas[j].y = 2*size_y - particulas[j].y;
                else if(particulas[j].y < 0)
                    particulas[j].y *= -1;
            }
        }
    }
    printf("\n\nSaida:\n%d %f %f\n", n_partic, size_x, size_y);
    for(i = 0; i < n_partic; ++i)
        printf("\t%f %f %f %f %f %f\n", particulas[i].x, particulas[i].y, particulas[i].vx, particulas[i].vy, particulas[i].carga, particulas[i].raio);
    return 0;
}
