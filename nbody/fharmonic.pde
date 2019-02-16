
void fharmonic(Particle[] ps) {
  for(Particle p : ps) {
    for(int k=0; k<3;++k) {
      p.hforce[k] += p.charge * ( -2.0 * voltage / rNullSquared ) * (p.pos[k] - rmin);
    }
  } 
}
