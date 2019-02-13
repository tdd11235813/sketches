import queasycam.*;

QueasyCam cam;

int np = 10;

double ts = 1e-2;//1e-14; // timestep in [s]
double particleMass = 4.03594014E-27+2*9.1093826E-31;// M(Mg)/A =  4.03594014⋅10−27 kg + 2me
// Mg(2+) Ion => 2e = 2x 1.6021766208(98)×10−19 C = 3.204353e-19 [C]
double particleCharge = 3.204353e-19; // [C] // 2*e
double voltage = 5; // [V]
double rmin    = 0; // position of potential minimum x=y=z=0
double rNull   = 1e-6; // [m]
double rNullSquared      = rNull*rNull; // [m]
double _md_phys_c        = 299792458.;
double _md_phys_emfactor = _md_phys_c*_md_phys_c*1E-7;
double EPS2    = 0.01;

class Particle {
  double[] pos;
  double[] vel = {0.0,0.0,0.0};
  double[] hforce = {0.0,0.0,0.0};
  double[] cforce = {0.0,0.0,0.0};
  double[] force = {0.0,0.0,0.0}; // total force
  double[] acc = {0.0,0.0,0.0}; // acceleration
  double mass;
  Particle() {
    pos = new double[3];
    pos[0] = randomGaussian();
    pos[1] = randomGaussian();
    pos[2] = randomGaussian();
    mass = particleMass;
  }
  void calc_force() {
    for(int k=0; k<3; ++k)
      force[k] = hforce[k] + cforce[k];
  }
  void calc_acc() {
    for(int k=0; k<3; ++k)
      acc[k] = force[k] / mass;
  }
  void integrate() {
    for(int k=0; k<3; ++k) {
      vel[k] += ts*acc[k];
      pos[k] += ts*vel[k];
    }
  }
};

Particle[] ps = new Particle[np];

void setup(){
  size(400, 400, P3D);
  
  cam = new QueasyCam(this);
  cam.speed = 2;              // default is 3
  cam.sensitivity = 0.25;      // default is 2
  println("Particles "+ps.length);
  for (int i = 0; i < ps.length; i++) {
    ps[i] = new Particle();
  }
}

void pcoulomb() {
  for(Particle p : ps) { 
    for(Particle q : ps) {
      double[] diff = { p.pos[0] - q.pos[0], p.pos[1] - q.pos[1], p.pos[2] - q.pos[2] };
      double distSqr = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
      double dist = java.lang.Math.sqrt(distSqr);
      double distCube = dist*distSqr;
      if(distCube>0.0) {
        double force_factor = _md_phys_emfactor * particleCharge * particleCharge / distCube;
        
        p.cforce[0] += force_factor * diff[0];
        p.cforce[1] += force_factor * diff[1];
        p.cforce[2] += force_factor * diff[2];
        q.cforce[0] += force_factor * diff[0];
        q.cforce[1] += force_factor * diff[1];
        q.cforce[2] += force_factor * diff[2]; 
      }
    }
  }
}

void pmove() {
  for(Particle p : ps) {
    p.calc_force();
    p.calc_acc();
    p.integrate();
  }
}

void draw(){
  pcoulomb();
  pmove();

  background(0);
  box(200);
  strokeWeight(1);
  scale(10);
  beginShape(POINTS);
  for(Particle p : ps) {
    vertex((float)p.pos[0],(float)p.pos[1],(float)p.pos[2]);
  }
  endShape();
}
