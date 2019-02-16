void fcoulomb(Particle[] ps) {
  for(Particle p : ps) {
    for(Particle q : ps) {
      if(p!=q) {
        double[] diff = { p.pos[0] - q.pos[0], p.pos[1] - q.pos[1], p.pos[2] - q.pos[2] };
        double distSqr = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
        double dist = java.lang.Math.sqrt(distSqr);
        double distCube = dist*distSqr;
        if(distCube>0.0) {
          double force_factor = _md_phys_emfactor * p.charge * p.charge / distCube;
          for(int k=0; k<3; ++k) {
            p.cforce[k] += force_factor * diff[k];
            q.cforce[k] -= force_factor * diff[k];
          }
        }
      }
    }
  }
}
