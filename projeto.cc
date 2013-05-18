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
    int i, j, k, tid, rodada = 0;
    std::ifstream file(argv[1], std::ifstream::in);
    if(!file) return 1;
    int n_partic = 0;
    double size_x = 0, size_y = 0;
    double c = atof(argv[2]);
    double epsilon = atof(argv[3]);
    double tau = atof(argv[4]);
    int iters = atoi(argv[5]);
    int procs = atoi(argv[6]);
    double fx, fy;
    
    omp_set_num_threads(procs);

    file >> n_partic;
    file >> size_x;
    file >> size_y;
    particula particulas[2][n_partic];
    Vector2D forcas[n_partic];
    int iterador[n_partic];
    for(i = 0; i < n_partic; ++i) {
        file >> particulas[rodada][i].x; file >> particulas[rodada][i].y;
        file >> particulas[rodada][i].vx; file >> particulas[rodada][i].vy;
        file >> particulas[rodada][i].carga; file >> particulas[rodada][i].raio;
    }
    printf("%s %f %f %f %d %d\n", argv[1], c, epsilon, tau, iters, procs);
    printf("%d, %f, %f\n", n_partic, size_x, size_y);
    for(i = 0; i < n_partic; ++i)
        printf("\t%f %f %f %f %f %f\n", particulas[rodada][i].x, particulas[rodada][i].y, particulas[rodada][i].vx, particulas[rodada][i].vy, particulas[rodada][i].carga,particulas[rodada][i].raio);
    printf("\n\n");
    
    #pragma omp parallel private(i, j, k, x, y, distancia, distancia_sqr, tid, rodada, fx, fy)
    {
        tid = omp_get_thread_num();
        x = 0; y = 0; fx = 0; fy = 0; rodada = 0; distancia = 0; distancia_sqr = 0; k = 0;
        for(i = 0; i < 10; ++i) {
            rodada = !rodada;
            //printf("%d",tid);
            //if(i%25000 == 0 || i > 2499900)
            //    printf("Iteracao %d\n", i);
            #pragma omp for
            for(j = 0; j < n_partic; ++j) {
                forcas[j].x(0);
                forcas[j].y(0);
                #pragma omp critical
                particulas[rodada][j].vx = particulas[!rodada][j].vx;
                particulas[rodada][j].vy = particulas[!rodada][j].vy;
                particulas[rodada][j].x = particulas[!rodada][j].x;
                particulas[rodada][j].y = particulas[!rodada][j].y;
                
                for(k = 0; k < n_partic; ++k) {
                    if(j == k)
                        continue;
                    x = (particulas[!rodada][k].x - particulas[rodada][j].x);
                    y = (particulas[!rodada][k].y - particulas[rodada][j].y);
                    distancia_sqr = x*x + y*y;
                    distancia = sqrt(distancia_sqr);
                    if(distancia > epsilon || distancia < particulas[!rodada][k].raio+particulas[0][j].raio)
                        continue;
                    forca = c/distancia_sqr;
                    fx = -(x/distancia_sqr)*forca*tau*particulas[0][j].carga*particulas[0][k].carga;
                    fy = -(y/distancia_sqr)*forca*tau*particulas[0][k].carga*particulas[0][k].carga;
                    forcas[j].add(fx, fy);
                }
                particulas[rodada][j].vx = particulas[!rodada][j].vx + forcas[j].x();
                particulas[rodada][j].vy = particulas[!rodada][j].vy + forcas[j].y();
            }
            
            /*if(tid == 0) {
                printf("Iter %d:\n", i);
                for(j = 0; j < n_partic; ++j)
                    printf("\t%f %f %f %f %f %f\n", particulas[rodada][j].x, particulas[rodada][j].y, particulas[rodada][j].vx, particulas[rodada][j].vy, particulas[rodada][j].carga, particulas[rodada][j].raio);
            }
            #pragma omp barrier
            */
            #pragma omp for
            for(j = 0; j < n_partic; ++j) {
                particulas[rodada][j].x += particulas[rodada][j].vx*tau;
                particulas[rodada][j].y += particulas[rodada][j].vy*tau;
                if(particulas[rodada][j].x > size_x) {
                    particulas[rodada][j].x = 2*size_x - particulas[rodada][j].x;
                    particulas[rodada][j].vx /= 2;
                }
                else if(particulas[rodada][j].x < 0) {
                    particulas[rodada][j].x *= -1;
                    particulas[rodada][j].vx /= 2;
                }
                
                if(particulas[rodada][j].y > size_y) {
                    particulas[rodada][j].y = 2*size_y - particulas[rodada][j].y;
                    particulas[rodada][j].vy /= 2;
                }
                else if(particulas[rodada][j].y < 0) {
                    particulas[rodada][j].y *= -1;
                    particulas[rodada][j].vy /= 2;
                }
                
            }
        }
    }
    printf("\n\nSaida:\n%d %f %f\n", n_partic, size_x, size_y);
    for(i = 0; i < n_partic; ++i)
        printf("\t%f %f %f %f %f %f\n", particulas[rodada][i].x, particulas[rodada][i].y, particulas[rodada][i].vx, particulas[rodada][i].vy, particulas[rodada][i].carga, particulas[rodada][i].raio);
    return 0;
}
