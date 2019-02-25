
class Particle {
  double[] pos;
  double[] vel = {0.0,0.0,0.0};
  double[] hforce = {0.0,0.0,0.0};
  double[] cforce = {0.0,0.0,0.0};
  double[] lforce = {0.0,0.0,0.0};
  double[] force = {0.0,0.0,0.0}; // total force
  double[] acc = {0.0,0.0,0.0}; // acceleration
  double mass;
  double charge;
  
  float vhforce;
  float vcforce;
  float vlforce;
  float vforce;
  PVector[] trajs;
  int ltrajs = 0; //current length of trajs
  
  Particle(double particleMass, double particleCharge, double temp) {
    pos = new double[3];
    trajs = new PVector[ntrajs];
    for(int s=0; s<ntrajs; ++s)
      trajs[s] = new PVector(0,0,0);

    pos[0] = normal_dist(0,ion_pos_sd);
    pos[1] = normal_dist(0,ion_pos_sd);
    pos[2] = normal_dist(0,ion_pos_sd);
    vel[0] = random_vel_temp(particleMass, temp);
    vel[1] = random_vel_temp(particleMass, temp);
    vel[2] = random_vel_temp(particleMass, temp);
    mass = particleMass;
    charge = particleCharge;
  }
  
  void calc_force_and_acc() {
    for(int k=0; k<3; ++k) {
      force[k] = hforce[k] + cforce[k] + lforce[k];
      acc[k]   = force[k] / mass;
    }
  }
  
  void integrate_euler(double ts) {
    for(int k=0; k<3; ++k) {
      vel[k] += ts*acc[k];
      pos[k] += ts*vel[k];
    }
  }
  
  void integrate_verlet_part1(double ts) {
    for(int k=0; k<3; ++k) {
      pos[k] += ts*(vel[k]+0.5*ts*acc[k]); // pos+=dt(v+a*dt2)
      vel[k] += 0.5*ts*acc[k]; // v += a*dt2
    }
  }
  void integrate_verlet_part2(double ts) {
    for(int k=0; k<3; ++k) {
      vel[k] += 0.5*ts*acc[k]; // v += a*dt2
    }
  }

  double get_v() {
    return dlen(vel);
  }
  double get_ekin() {
    return 0.5*mass*get_v()*get_v(); // [J]
  }

  void update_trajectory() {
    PVector npos = new PVector((float)pos[0],(float)pos[1],(float)pos[2]);
    if(ltrajs>0 && npos.dist(trajs[ltrajs-1])<5e-6)
      return; // do not add because distance to last point to small
    ++ltrajs;
    if(ltrajs>ntrajs) {
      ltrajs=ntrajs;
      for(int s=1;s<ntrajs;++s) {
        trajs[s-1] = trajs[s]; 
      }
    }
    trajs[ltrajs-1] = npos;
  }
  void reset_forces() {
    for(int k=0; k<3; ++k) {
      cforce[k] = 0;
      hforce[k] = 0;
      lforce[k] = 0;
    }
  }
};
