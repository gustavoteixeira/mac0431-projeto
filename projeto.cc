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
    int i, j, k, tid, rodada = 0;
    std::ifstream file(argv[1], std::ifstream::in);
    std::ofstream output("saida.txt");
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
    particula particulas[n_partic];
    Vector2D forcas[n_partic];
    int iterador[n_partic];
    for(i = 0; i < n_partic; ++i) {
        file >> particulas[i].x; file >> particulas[i].y;
        file >> particulas[i].lx; file >> particulas[i].ly;
        file >> particulas[i].carga; file >> particulas[i].raio;
    }
    
    #pragma omp parallel private(i, j, k, x, y, distancia, distancia_sqr, tid, rodada, fx, fy)
    {
        //tid = omp_get_thread_num();
        x = 0; y = 0; fx = 0; fy = 0; rodada = 0; distancia = 0; distancia_sqr = 0; k = 0;
        for(i = 0; i < iters; ++i) {
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
                    if(distancia > epsilon || distancia < particulas[k].raio + particulas[j].raio)
                        continue;
                    forca = c/distancia_sqr;
                    fx = -(x/distancia_sqr)*forca*tau*particulas[j].carga*particulas[k].carga;
                    fy = -(y/distancia_sqr)*forca*tau*particulas[j].carga*particulas[k].carga;
                    forcas[j].add(fx, fy);
                }
                particulas[j].lx = particulas[j].vx + forcas[j].x();
                particulas[j].ly = particulas[j].vy + forcas[j].y();
            }
            
            #pragma omp for
            for(j = 0; j < n_partic; ++j) {
                particulas[j].x += particulas[j].vx*tau;
                particulas[j].y += particulas[j].vy*tau;
                if(particulas[j].x > size_x) {
                    particulas[j].x = 1.5*size_x - 0.5*particulas[j].x;
                    particulas[j].vx *= -0.5;
                }
                else if(particulas[j].x < 0) {
                    particulas[j].x *= -0.5;
                    particulas[j].vx *= -0.5;
                }
                
                if(particulas[j].y > size_y) {
                    particulas[j].y = 1.5*size_y - 0.5*particulas[j].y;
                    particulas[j].vy *= -0.5;
                }
                else if(particulas[j].y < 0) {
                    particulas[j].y *= -0.5;
                    particulas[j].vy *= -0.5;
                }
                
            }
        }
    }
    
    output.precision(6);
    output << n_partic << " " << size_x << " " << size_y << std::endl;
    for(i = 0; i < n_partic; ++i)
        output << "    " << particulas[i].x  << " " << particulas[i].y  << " " <<
                            particulas[i].vx << " " << particulas[i].vy << " " <<
                            particulas[i].carga << " " << particulas[i].raio << std::endl;
        
    output.close();
    return 0;
}
