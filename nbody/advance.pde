
void advance(Particle[] ps, double ts) {
  for(Particle p : ps) {
    p.calc_force();
    p.calc_acc();
    p.integrate(ts);
  }
}
