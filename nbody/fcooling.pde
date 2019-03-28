
double cooling_linear(double vk, double ts, double mass) {
    double p=0.1; // 10%
    double dv = vk;
    double ln1p = java.lang.Math.log(1+p)/java.lang.Math.log(java.lang.Math.exp(1.0));
    double restore = -ln1p/ts * mass;  
    return restore * dv;
}

void fcooling_linear(Particle[] ps, double ts) {
  for(Particle p : ps) {
    for(int k=0; k<3; ++k)
      p.lforce[k] += cooling_linear(p.vel[k], ts, p.mass);
  }
}

void fcooling_russian(Particle[] ps) {
  double vmax = 0.0;
  for(Particle p : ps) { 
    double vm = p.get_v();
    //if(vm>vmax)
    //  vm=vmax;
    if(vm>90) { // TODO:
      p.vel[0] = 0;
      p.vel[1] = 0;
      p.vel[2] = 0;
    }
  }
}
