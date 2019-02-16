
double restore = 1e-16; // [ C*V/m*s/m ], too soft, v still too high 

double cooling_linear(double vk, double vacc, double vmax) {
    double dv = vk - vmax;
    if(vacc <0.) return -restore * dv;
    else if ((dv < 0.) && (dv > -vacc)) return +restore * (dv+vacc);
    else if ((dv > 0.) && (dv < +vacc)) return -restore * (dv-vacc);
    else return 0.0;
}

// F = 0                   <=>        |v-vmax| > vacc                    //
// F = - D * (v-vmax-vacc) <=>     0 < v-vmax < vacc                     //
// F = + D * (v-vmax+vacc) <=> -vacc < v-vmax < 0                        //
void fcooling_linear(Particle[] ps, double vacc, double vmax) {
  for(Particle p : ps) {
    for(int k=0; k<3; ++k)
      p.lforce[k] += cooling_linear(p.vel[k], vacc, vmax);
  }
}

void fcooling_russian(Particle[] ps) {
  double vmax = 0.0;
  for(Particle p : ps) { 
    double vm = java.lang.Math.sqrt(p.vel[0]*p.vel[0]+p.vel[1]*p.vel[1]+p.vel[2]*p.vel[2]);
    //if(vm>vmax)
    //  vm=vmax;
    if(vm>1e+3) { // TODO:
      p.vel[0] = 0;
      p.vel[1] = 0;
      p.vel[2] = 0;
    }
  }
}