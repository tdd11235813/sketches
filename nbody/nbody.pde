import queasycam.*;

// TODO
// - init temp -> v

boolean show_dist_markers = false;
boolean show_trajectories = false;
boolean use_velocity_verlet = true; // velocity-verlet integrator
double init_vel_by_random_temp = 10; // [K] // initial velocity by random temperature distribution

int np = 500;
int ntrajs = 100;

double ts                = 5e-8;//1e-14; // timestep in [s]
double particleMass      = 4.03594014E-27+2*9.1093826E-31;// M(Mg)/A =  4.03594014⋅10−27 kg + 2me
// Mg(2+) Ion => 2e = 2x 1.6021766208(98)×10−19 C = 3.204353e-19 [C]
double particleCharge    = 3.204353e-19; // [C] // 2*e
double voltage           = 5; // [V]
double rmin              = 0; // position of potential minimum x=y=z=0
double rNull             = 2e-2; // [m]
double rNullSquared      = rNull*rNull; // [m]
double _md_phys_c        = 299792458.;
double _md_phys_emfactor = _md_phys_c*_md_phys_c*1E-7;
double _md_phys_mu0      = PI*4.E-7;
double _kb               = 1.38064852e-23; // [J K^-1], [m^2⋅kg/(s^2⋅K)]
double eps0              = 1./(_md_phys_mu0*_md_phys_c*_md_phys_c);
double ion_pos_sd        = 5e-4; // [m]

double ctime = 0;

Particle[] ps = new Particle[np];

QueasyCam cam;
PMatrix3D baseMat;
// ---------------------------------------------------------

void setup(){
  //size(400, 400, P3D);
  fullScreen(P3D);
  noCursor();
  hint(DISABLE_DEPTH_SORT);
  hint(DISABLE_DEPTH_TEST);
  // Remember the start model view matrix values
  baseMat = getMatrix(baseMat);
  
  cam = new QueasyCam(this);
  cam.speed = 0.1;              // default is 3
  cam.sensitivity = 0.1;      // default is 2
  cam.position = new PVector(20,-20,-40);
  cam.pan = 2*PI/3;
  cam.tilt = PI/8;
  println("Particles "+ps.length);
  for (int i = 0; i < ps.length; i++) {
    ps[i] = new Particle(particleMass, particleCharge, init_vel_by_random_temp);
  }
}

// ---------------------------------

PVector gridPos = new PVector(0, 0, 0);//the position of the grid (it moves with the camera to make it look infinite)

void rectGrid(int size, int tilesize, float y) {
  noFill();//i only want the outline of the rectangles
  for (float x = -size/2; x <= size/2; x++) {
    for (float z = -size/2; z <= size/2; z++) {
      //run two for loops, cycling through 10 different positions of rectangles
      pushMatrix();

      stroke(200, 200, 200, 0.25*map(dist(-gridPos.x, -gridPos.z, x*tilesize, z*tilesize), 0, size/2*tilesize, 255, 0));//the rectangles close to you, are clear, while the ones farther from you, are much fainter
      //uncomment the next line:
      //stroke(0,255,0);
      // to see how the infinity thing works

      translate(x*tilesize, y, z*tilesize);//move the rectangles to where they shall be
      rotateX(HALF_PI);
      rect(0, 0, tilesize, tilesize);
      popMatrix();
    }
  }
}

// ---------------------------------

void draw(){
  // scale positions up
  double sc = 2e4;
  double dist_thr = 1.1e-4;
  float sc_force = 25e+5;
  // -- forces --
  // init coulomb force to 0
  for(Particle p : ps) {
    p.reset_forces();
  }
  fcoulomb(ps);
  fharmonic(ps);
  //fcooling_linear(ps,-1,1);
  fcooling_linear(ps,-20,-10);
  fcooling_linear(ps, 20, 10);
  //fcooling_russian(ps);
  // -- integrate and move --
  // advance(ps,ts); // euler integrator
  advance_verlet(ps,ts); // velocity-verlet integrator
  ctime += ts;
  // -- render scene --

  background(250);
  strokeWeight(1);
  rectGrid(40,2,0);

  if(show_dist_markers) {
    // distance marker
    strokeWeight(1.5);
    beginShape(LINES);
    for(Particle p : ps) {
      for(Particle q : ps) {
        double[] diff = { p.pos[0] - q.pos[0], p.pos[1] - q.pos[1], p.pos[2] - q.pos[2] };
        double distSqr = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
        double dist = java.lang.Math.sqrt(distSqr);
        if(dist < dist_thr) {
          stroke(0, (float)((dist_thr-dist)*125/dist_thr));
          PVector pos1=new PVector(
          (float)(sc*p.pos[0]),
          (float)(sc*p.pos[1]),
          (float)(sc*p.pos[2]));
          PVector pos2=new PVector(
          (float)(sc*q.pos[0]),
          (float)(sc*q.pos[1]),
          (float)(sc*q.pos[2]));
          vertex(pos1.x,pos1.y,pos1.z);
          vertex(pos2.x,pos2.y,pos2.z);
        }
      }
    }
    endShape();
  }
  
  if(show_trajectories) {
    pushMatrix();
    scale((float)sc);
    strokeWeight(1.0/(float)sc);
    beginShape(LINES);
    for(Particle p : ps) {
      p.update_trajectory();
      for(int s=1; s<p.ltrajs; ++s) {
        PVector x0 = p.trajs[s-1];
        PVector x1 = p.trajs[s];
        stroke(0,230.0*s/ntrajs+25);
        vertex(x0.x,x0.y,x0.z);
        vertex(x1.x,x1.y,x1.z);
      }
    }
    endShape();
    popMatrix();
  }
  // render particles as points
  strokeWeight(11);
  beginShape(POINTS);
  for(Particle p : ps) {
    PVector pos=new PVector(
        (float)(sc*p.pos[0]),
        (float)(sc*p.pos[1]),
        (float)(sc*p.pos[2]));
    float z = PVector.dist(pos,cam.position);
    stroke(
      0,0,0, 
      max(255-220.0*pow(z/50,4),25));
    vertex(pos.x,pos.y,pos.z);
  }
  endShape();
  strokeWeight(8);
  beginShape(POINTS);
  for(Particle p : ps) {
    PVector pos=new PVector(
        (float)(sc*p.pos[0]),
        (float)(sc*p.pos[1]),
        (float)(sc*p.pos[2]));
    float z = PVector.dist(pos,cam.position);
    stroke(
      sc_force*p.vcforce/np, 
      sc_force*p.vhforce/np, 
      sc_force*p.vlforce/np, 
      max(255-220.0*pow(z/50,4),25));
    vertex(pos.x,pos.y,pos.z);
  }
  endShape();
  
  this.setMatrix(baseMat);
  ambientLight(255, 255, 255); 
  show_gui();
}
