
void advance(Particle[] ps, double ts) {
  for(Particle p : ps) {
    p.calc_force();
    p.calc_acc();
    p.integrate(ts);
  }
}

void advance_verlet(Particle[] ps, double ts) {
  for(Particle p : ps) {
    p.calc_force();
    p.calc_acc();
    p.integrate_verlet1(ts);
    p.calc_force();
    p.calc_acc();
    p.integrate_verlet2(ts);
  }
}
