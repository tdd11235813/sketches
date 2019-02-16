

double dlen(double[] vec3d) {
  double distSqr = vec3d[0]*vec3d[0] + vec3d[1]*vec3d[1] + vec3d[2]*vec3d[2];
  return java.lang.Math.sqrt(distSqr);
}
double ddist(double[] vec3d1, double[] vec3d2) {
    double[] diff = { vec3d1[0] - vec3d2[0], vec3d1[1] - vec3d2[1], vec3d1[2] - vec3d2[2] };
    return dlen(diff);
}

// ---

float avg_min_distance() {
  double avg_min_dist = 0;
  for( Particle p : ps ) {
    double min_dist = 1e9;
    for( Particle q : ps ) {
      if(p!=q) {
        double dist = ddist(p.pos, q.pos);
        if( dist < min_dist )
          min_dist = dist;
      }
    }
    if( min_dist < 1e9 )
      avg_min_dist += min_dist; 
  }
  return (float)(avg_min_dist / ps.length);
}

float avg_v() {
  double avg_v = 0;
  for( Particle p : ps ) {
    avg_v += p.get_v();
  }
  return (float)(avg_v/ps.length);
}

float avg_ekin() {
  double avg_ekin = 0;
  for( Particle p : ps ) {
    avg_ekin += p.get_ekin();
  }
  return (float)(avg_ekin/ps.length);
}

float avg_temp() {
  double avg_temp = 0;
  for( Particle p : ps ) {
    double temp = 2*p.get_ekin()/_kb;
    avg_temp += temp;
  }
  return (float)(avg_temp / ps.length);
}

double get_ekin() {
  double ekin = 0;
  for( Particle p : ps ) {
    ekin += p.get_ekin();
  }
  return ekin;
}

double get_epot() {
  double epot=0;
  for( Particle p : ps ) {
    
    // epot in harmonic potential
    double dp = dlen(p.pos); // ||position - minpos||
    double epoth = p.charge * voltage / rNullSquared * dp * dp; // [C*V/m² * m] = [J/m] = [N]

    // epot in coulomb potential
    double epotc = 0;
    for( Particle q : ps ) {
      if(p!=q) {
        epotc += 1/(4*PI*eps0) * p.charge*q.charge/(ddist(p.pos,q.pos));
      }
    }
    
    epot += epotc + epoth;
  }
  return epot;
}
