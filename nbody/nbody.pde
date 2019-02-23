import queasycam.*;

boolean show_dist_markers = false;
boolean show_trajectories = false;
boolean use_velocity_verlet = true; // velocity-verlet integrator
double init_vel_by_random_temp = 10; // [K] // initial velocity by random temperature distribution

int np                   = 500; // number of ions
int ntrajs               = 100; // trajectory steps
double ion_pos_sd        = 5e-4; // [m] initial position range, sd for gaussian distribution

double ts                = 1e-7; // timestep in [s]
// using ^24 Mg^+ ions
double particleMass      = 24*1.66053886E-27; // [C]
double particleCharge    = 1.6021766209e-19;  // [C]
double voltage           = 5; // [V]
double rmin              = 0; // position of potential minimum x=y=z=0
double rNull             = 1.5e-2; // [m]
double rNullSquared      = rNull*rNull; // [m]

double _md_phys_c        = 299792458.;
double _md_phys_emfactor = _md_phys_c*_md_phys_c*1E-7;
double _md_phys_mu0      = PI*4.E-7;
double _kb               = 1.38064852e-23; // [J K^-1], [m^2⋅kg/(s^2⋅K)]
double _md_phys_h        = 6.6260693E-34;
double eps0              = 1./(_md_phys_mu0*_md_phys_c*_md_phys_c);

double ctime = 0;

Particle[] ps = new Particle[np];

// currently not used, because it is not working
CoolingLaser laser1      = new CoolingLaser(280.0,  0.57735,0.57735,0.57735); // wavelength [nm], unit direction vector
CoolingLaser laser2      = new CoolingLaser(280.0, -0.57735,-0.57735,-0.57735);

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
  double sc = 5e4;
  double dist_thr = 6e-5;
  float sc_force = 50e+5;
  // -- forces --
  // init forces to 0
  for(Particle p : ps) {
    p.reset_forces();
  }
  // Coulomb interaction
  fcoulomb(ps);
  // harmonic oscillator
  fharmonic(ps);
  // cooling
  //fcooling_linear(ps,-1, 1);
  //fcooling_linear(ps,-1,-1);
  fcooling_linear(ps,-20,-10);
  fcooling_linear(ps, 20, 10);
  //fcooling_russian(ps);
  //fcooling(ps, laser1);
  //fcooling(ps, laser2);
  
  // -- integrate and move --
  if(use_velocity_verlet) {
    advance_verlet(ps,ts); // velocity-verlet integrator
  } else {
    advance(ps,ts); // euler integrator
  }

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
    strokeWeight(1.0);
    beginShape(LINES);
    for(Particle p : ps) {
      p.update_trajectory();
      for(int s=1; s<p.ltrajs; ++s) {
        PVector x0 = PVector.mult(p.trajs[s-1],(float)sc);
        PVector x1 = PVector.mult(p.trajs[s],(float)sc);
        float z = PVector.dist(x1,cam.position);
        stroke(0,(230.0*s/ntrajs+25)*max(1.0-0.9*pow(z/50,4),0.1));
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
