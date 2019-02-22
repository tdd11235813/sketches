// https://pdfs.semanticscholar.org/3d80/7a3327c1768a030dfc96345f6f6263358c13.pdf
class CoolingLaser {
  double wavelength; // [m]
  double G = TWO_PI*42.7e6; // natural linewidth of the force profile [Hz]
  double S = 3; // saturation parameter
  double detuning = 6.5*G; // [Hz] // wL - w0

  double dir[] = {0.,0.,0.}; // direction (k), unit vector
  
  CoolingLaser( double wavelength_nm, double kx, double ky, double kz ) {
    wavelength = wavelength_nm*1e-9;
    dir[0] = kx;
    dir[1] = ky;
    dir[2] = kz;
  }
  // |k| = |dir*2Pi/wavelength| = 2Pi/wavelength
  // F_L = hbar/8 * ( |k|*S*G^3/( (detuning-dot(v,k)) )^2+G^2/4*(1+S) )
  void cooling(Particle[] ps) {
    double Gdiv2     = G / 2.;
    double Gdiv2sqr  = Gdiv2 * Gdiv2;
    double Gdiv2cube = Gdiv2sqr * Gdiv2;
    double t1  = _md_phys_h * S * Gdiv2cube / wavelength;    
    double t2  = Gdiv2sqr * ( 1 + S );
    
    for(Particle p : ps) {
      double vk = 0;
      for(int k=0; k<3; ++k) {
        vk += dir[k] * p.vel[k];
      }
      double t3 = (detuning - vk*TWO_PI/wavelength);
      t3 *= t3;
      double t4 = t1 / (t3+t2);
      for(int k=0; k<3; ++k) {        
        p.lforce[k] += t4 * dir[k]; // +?
      }
    }
  }
};
