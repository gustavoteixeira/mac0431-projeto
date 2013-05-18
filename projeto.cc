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
    double x, y, vx, vy, lx, ly, carga, raio;
};

typedef struct part particula;

int main(int argc, char* argv[]) {
    double x, y, distancia_sqr, distancia, forca = 0;
    int i, j, k, tid;
    bool bounce = false;
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
    for(i = 0; i < n_partic; ++i) {
        file >> particulas[i].x; file >> particulas[i].y;
        file >> particulas[i].lx; file >> particulas[i].ly;
        file >> particulas[i].carga; file >> particulas[i].raio;
    }
    printf("%s %f %f %f %d %d\n", argv[1], c, epsilon, tau, iters, procs);
    printf("%d, %f, %f\n", n_partic, size_x, size_y);
    for(i = 0; i < n_partic; ++i)
        printf("\t%f %f %f %f %f %f\n", particulas[i].x, particulas[i].y, particulas[i].vx, particulas[i].vy, particulas[i].carga,particulas[i].raio);
    printf("\n\n");
    
    #pragma omp parallel private(i, j, k, x, y, distancia, distancia_sqr, bounce, tid)
    {
        tid = omp_get_thread_num();
        x = 0; y = 0; distancia = 0; distancia_sqr = 0; k = 0; bounce = false;
        for(i = 0; i < iters; ++i) {
            //printf("%d",tid);
            //if(i%25000 == 0 || i > 2499900)
            //    printf("Iteracao %d\n", i);
            #pragma omp for
            for(j = 0; j < n_partic; ++j) {
                forcas[j].x(0);
                forcas[j].y(0);
                particulas[j].vx = particulas[j].lx;
                particulas[j].vy = particulas[j].ly;
                for(k = 0; k < n_partic; ++k) {
                    if(j == k)
                        continue;
                    x = (particulas[k].x - particulas[j].x);
                    y = (particulas[k].y - particulas[j].y);
                    distancia_sqr = x*x + y*y;
                    distancia = sqrt(distancia_sqr);
                    if(distancia > epsilon)
                        continue;
                    if(distancia < particulas[k].raio+particulas[j].raio)
                        continue;
                    forca = c/distancia_sqr;
                    forcas[j].add(-(x/distancia_sqr)*forca*tau*particulas[j].carga*particulas[k].carga, -(y/distancia_sqr)*forca*tau*particulas[k].carga*particulas[k].carga);
                }
                particulas[j].lx = particulas[j].vx + forcas[j].x();
                particulas[j].ly = particulas[j].vy + forcas[j].y();
            }
            
            /*if(tid == 0) {
                printf("Iter %d:\n", i);
                for(j = 0; j < n_partic; ++j)
                    printf("\t%f %f %f %f %f %f\n", particulas[j].x, particulas[j].y, particulas[j].vx, particulas[j].vy, particulas[j].carga, particulas[j].raio);
            }*/
            //#pragma omp barrier
            #pragma omp for
            for(j = 0; j < n_partic; ++j) {
                particulas[j].x += particulas[j].vx*tau;
                particulas[j].y += particulas[j].vy*tau;
                if(particulas[j].x > size_x) {
                    particulas[j].x = 2*size_x - particulas[j].x;
                    bounce = true;
                }
                else if(particulas[j].x < 0) {
                    particulas[j].x *= -1;
                    bounce = true;
                }
                
                if(particulas[j].y > size_y) {
                    particulas[j].y = 2*size_y - particulas[j].y;
                    bounce = true;
                }
                else if(particulas[j].y < 0) {
                    particulas[j].y *= -1;
                    bounce = true;
                }
                if(bounce) {
                    bounce = false;
                    particulas[j].vx /= 2;
                    particulas[j].vy /= 2;
                }
            }
        }
    }
    printf("\n\nSaida:\n%d %f %f\n", n_partic, size_x, size_y);
    for(i = 0; i < n_partic; ++i)
        printf("\t%f %f %f %f %f %f\n", particulas[i].x, particulas[i].y, particulas[i].vx, particulas[i].vy, particulas[i].carga, particulas[i].raio);
    return 0;
}
