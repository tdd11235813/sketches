
float amdist = 0;
double ekin = 0;
double epot = 0;
double temp = 0;
double vavg_v = 0;
double vavg_temp = 0;
double vavg_ekin = 0;
double gamma = 0;

int frame = 0;

void show_gui() {
  if((frame%10)==0) {
    ekin = get_ekin();
    epot = get_epot();
    amdist = avg_min_distance()*1e6; // [m] -> [µm]
    vavg_v = avg_v();
    vavg_ekin = avg_ekin();
    vavg_temp = avg_temp();
    temp = 2.0/3.0*ekin/np*1.0/_kb;
    gamma = epot/ekin;
    frame = 0;
  }
  ++frame;
  
  textAlign(RIGHT);
  textSize(15);

  int toffset = 175;
  int voffset = 55;
  int yoffset = 25;
  pushMatrix();
  stroke(125);
  strokeWeight(1);
  fill(233,244,255,175);
  translate(width-100,100);
  rect(-330,-35,375, 375);
  noStroke();
  fill(50,50,100);

  text("Ions (Mg2+): ", -toffset, 0); 
   text(ps.length, -voffset, 0);
  translate(0,yoffset);
  text("Mass: ", -toffset, 0); 
   text(String.format("%g",particleMass), -voffset, 0); text("[kg]", 0, 0);
  translate(0,yoffset);
  text("Charge: ", -toffset, 0); 
   text(String.format("%g",particleCharge), -voffset, 0); text("[C]", 0, 0);
  translate(0,yoffset);
  text("TimeStep: ", -toffset, 0); 
   text(String.format("%g",ts), -voffset, 0); text("[s]", 0, 0);
  translate(0,yoffset);
  text("Time: ", -toffset, 0); 
   text(String.format("%g",ctime*1e6), -voffset, 0); text("[µs]", 0, 0);

  translate(0,1.5*yoffset);
  text("AvgMinDistance: ", -toffset, 0); text("[µm]", 0, 0);
   text(String.format("%g",amdist), -voffset, 0);
  translate(0,yoffset);
  text("AvgVel: ", -toffset, 0); text("[m/s]", 0, 0);
   text(String.format("%g",vavg_v), -voffset, 0);
  translate(0,yoffset);
  text("AvgEkin: ", -toffset, 0); text("[J]", 0, 0);
   text(String.format("%g",vavg_ekin), -voffset, 0);
  translate(0,yoffset);
  text("Ekin: ", -toffset, 0); text("[J]", 0, 0);
   text(String.format("%g",ekin), -voffset, 0);
  translate(0,yoffset);
  text("Epot: ", -toffset, 0); text("[J]", 0, 0);
   text(String.format("%g",epot), -voffset, 0);
  translate(0,yoffset);
  text("E: ", -toffset, 0); text("[J]", 0, 0);
   text(String.format("%g",epot+ekin), -voffset, 0);
  translate(0,yoffset);
  text("AvgTemp: ", -toffset, 0);
   text(String.format("%g",vavg_temp), -voffset, 0); text("[K]", 0, 0);
  translate(0,yoffset);
  text("GasTemp: ", -toffset, 0);
   text(String.format("%g",temp), -voffset, 0); text("[K]", 0, 0);
  translate(0,yoffset);
  text("Gamma: ", -toffset, 0);
   text(String.format("%g",gamma), -voffset, 0); text("[-]", 0, 0);
  popMatrix();
}
