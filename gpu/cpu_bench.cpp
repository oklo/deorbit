// CPU benchmark counterpart to gpu_run: same config, same fixed dt, K steps,
// reports ms/step + energy drift. (Plain C++; the existing parallel CPU solver.)
//   ./cpu_bench <L> <steps>
#include "../sph/sph.hpp"
#include <cstdio>
#include <chrono>
#include <random>

static double now(){ return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

int main(int argc,char** argv){
    int Lc = argc>1?atoi(argv[1]):60; int steps = argc>2?atoi(argv[2]):10;
    System S; S.materials.push_back(Material::basalt()); S.materials.push_back(Material::iron());
    S.eta=1.3; S.eps_as=0.3; S.damage_on=true;
    const double dx=30.0, R=0.2*Lc*dx; S.h_init=S.eta*dx;
    double lo0=-Lc*dx/2.0;
    S.set_domain(lo0-dx, -lo0+dx, false, lo0-dx, -lo0+dx, false, lo0-dx, -lo0+dx, false);
    std::mt19937 rng(1); std::uniform_real_distribution<double> J(-0.2*dx,0.2*dx);
    for(int a=0;a<Lc;a++)for(int b=0;b<Lc;b++)for(int c=0;c<Lc;c++){
        double x=lo0+(a+0.5)*dx+J(rng),y=lo0+(b+0.5)*dx+J(rng),z=lo0+(c+0.5)*dx+J(rng);
        bool iron=(x*x+y*y+z*z)<R*R; S.cur_mat=iron?1:0;
        S.add({x,y,z},{80.0*(2*b/(double)Lc-1),0.0,-40.0*(2*c/(double)Lc-1)},
              (iron?7800.0:2700.0)*dx*dx*dx, 2e5, false); }
    S.seed_damage(dx); S.init();
    auto energy=[&]{ double KE=0,IE=0; for(int i=0;i<S.n;i++){ KE+=0.5*S.mass[i]*norm2(S.vel[i]); IE+=S.mass[i]*S.u[i]; } return KE+IE; };
    double E0=energy(); double dt=2.0e-4;
    double t0=now(); for(int s=0;s<steps;s++) S.step(dt); double el=now()-t0;
    double vmax=0; for(int i=0;i<S.n;i++) vmax=std::max(vmax,norm(S.vel[i]));
    printf("CPU n=%d (%dx^3)  %d steps  %.1f ms/step  %.0f Mparticle-steps/s   E drift=%.2e  vmax=%.0f\n",
           S.n, Lc, steps, 1000*el/steps, S.n/(el/steps)/1e6, (energy()-E0)/E0, vmax);
    return 0;
}
