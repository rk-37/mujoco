// Compare MuJoCo's built-in ellipsoid fluid model against a real airfoil
// polar (with stall) injected via the mjcb_passive callback.
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mujoco/mujoco.h>

#define DEG (M_PI/180.0)

// ---- Airfoil polar: linear lift up to stall, then flat-plate post-stall ----
// alpha in radians. Returns lift & drag coefficients.
static void airfoil_polar(double alpha, double AR, double* Cl, double* Cd) {
  const double a_stall = 14*DEG;     // stall angle
  const double Cl_alpha = 2*M_PI;    // thin-airfoil lift slope (per rad)
  const double Cd0 = 0.02;           // parasitic drag
  const double e = 0.9;              // span efficiency

  // attached flow (linear lift, parabolic drag polar w/ induced drag)
  double Cl_lin = Cl_alpha * alpha;
  double Cd_lin = Cd0 + Cl_lin*Cl_lin/(M_PI*AR*e);
  // fully separated flow (flat plate)
  double Cl_sep = 2*sin(alpha)*cos(alpha);
  double Cd_sep = Cd0 + 2*sin(alpha)*sin(alpha);
  // smooth blend across stall
  double s = 1.0/(1.0 + exp(-(fabs(alpha) - a_stall)/(3*DEG)));
  *Cl = (1-s)*Cl_lin + s*Cl_sep;
  *Cd = (1-s)*Cd_lin + s*Cd_sep;
}

// ---- mjcb_passive callback: applies polar aero to named lifting surfaces ----
static int   aero_geom[8], n_aero = 0;
static double aero_S[8], aero_AR[8];

static void polar_passive(const mjModel* m, mjData* d) {
  double rho = m->opt.density;
  for (int i = 0; i < n_aero; i++) {
    int g = aero_geom[i];
    int body = m->geom_bodyid[g];
    // geom local axes (columns of geom_xmat): chord=x, span=y, normal=z
    const mjtNum* R = d->geom_xmat + 9*g;
    double cx[3] = {R[0], R[3], R[6]};   // chord axis
    double sy[3] = {R[1], R[4], R[7]};   // span axis
    double nz[3] = {R[2], R[5], R[8]};   // normal axis
    // geom velocity (world), relative to the medium (wind)
    mjtNum vel6[6];
    mj_objectVelocity(m, d, mjOBJ_GEOM, g, vel6, 0);
    double Vr[3] = { vel6[3]-m->opt.wind[0],
                     vel6[4]-m->opt.wind[1],
                     vel6[5]-m->opt.wind[2] };
    double V = sqrt(Vr[0]*Vr[0]+Vr[1]*Vr[1]+Vr[2]*Vr[2]);
    if (V < 1e-4) continue;

    double v_chord = Vr[0]*cx[0]+Vr[1]*cx[1]+Vr[2]*cx[2];
    double v_norm  = Vr[0]*nz[0]+Vr[1]*nz[1]+Vr[2]*nz[2];
    double alpha = atan2(-v_norm, v_chord);

    double Cl, Cd;
    airfoil_polar(alpha, aero_AR[i], &Cl, &Cd);
    double q = 0.5*rho*V*V, S = aero_S[i];
    double L = q*S*Cl, D = q*S*Cd;

    // drag opposes motion; lift perpendicular, in the chord-normal plane
    double drag_dir[3] = {-Vr[0]/V, -Vr[1]/V, -Vr[2]/V};
    double lift_dir[3];                       // = normalize(Vr x span)
    mju_cross(lift_dir, Vr, sy);
    double ln = sqrt(lift_dir[0]*lift_dir[0]+lift_dir[1]*lift_dir[1]+lift_dir[2]*lift_dir[2]);
    if (ln < 1e-6) continue;
    for (int k=0;k<3;k++) lift_dir[k] /= ln;

    mjtNum force[3] = { drag_dir[0]*D + lift_dir[0]*L,
                        drag_dir[1]*D + lift_dir[1]*L,
                        drag_dir[2]*D + lift_dir[2]*L };
    mjtNum torque[3] = {0,0,0};
    mj_applyFT(m, d, force, torque, d->geom_xpos+3*g, body, d->qfrc_passive);
  }
}

// glide a model to the ground, return horizontal distance & peak AoA seen
static double glide(mjModel* m, int use_polar, double* peak_aoa) {
  mjData* d = mj_makeData(m);
  mj_resetDataKeyframe(m, d, 0);
  mjcb_passive = use_polar ? polar_passive : NULL;
  *peak_aoa = 0;
  while (d->time < 12.0 && d->qpos[2] > 0.05) {
    // track AoA of the main wing (geom 0 in aero list) for reporting
    if (use_polar && n_aero) {
      mjtNum v6[6]; mj_objectVelocity(m,d,mjOBJ_GEOM,aero_geom[0],v6,0);
      const mjtNum* R=d->geom_xmat+9*aero_geom[0];
      double vc=v6[3]*R[0]+v6[4]*R[3]+v6[5]*R[6];
      double vn=v6[3]*R[2]+v6[4]*R[5]+v6[5]*R[8];
      double a=fabs(atan2(-vn,vc))/DEG; if(a>*peak_aoa)*peak_aoa=a;
    }
    mj_step(m, d);
  }
  double dist = d->qpos[0];
  mjcb_passive = NULL;
  mj_deleteData(d);
  return dist;
}

int main(void) {
  char err[1000]="";

  // ---------- 1) Static AoA sweep: ellipsoid vs airfoil polar ----------
  // single wing matching the glider's main wing: chord 0.24, span 1.4
  mjModel* mw = mj_loadXML("/tmp/wing.xml", NULL, err, sizeof(err));
  if(!mw){printf("load wing: %s\n",err);return 1;}
  double S  = M_PI*0.12*0.70;          // ellipse planform (chord/2 * span/2)
  double AR = 1.40/0.24;               // span/chord
  double Vinf = 15.0;
  printf("Wing: chord 0.24 m, span 1.40 m, area %.3f m^2, AR %.1f, V=15 m/s\n", S, AR);
  printf("\n        |------ ellipsoid (built-in) ------|------ airfoil polar (w/ stall) ------|\n");
  printf("  AoA     lift(N)  drag(N)   L/D    |   lift(N)  drag(N)   L/D     Cl\n");
  for (double aoa=0; aoa<=28.001; aoa+=2) {
    // ellipsoid via MuJoCo
    mjData* d = mj_makeData(mw);
    double a=aoa*DEG;
    d->qpos[3]=cos(a/2); d->qpos[5]=sin(a/2);
    mw->opt.wind[0]=Vinf; mju_zero(d->qvel,mw->nv);
    mj_forward(mw,d);
    double eL=d->qpos[2]>-1? d->qfrc_passive[2]:0, eD=d->qfrc_passive[0];
    mj_deleteData(d);
    // polar via formula
    double Cl,Cd; airfoil_polar(aoa*DEG, AR, &Cl,&Cd);
    double q=0.5*mw->opt.density*Vinf*Vinf;
    double pL=q*S*Cl, pD=q*S*Cd;
    printf("  %4.1f   %8.2f %7.2f %6.2f   | %8.2f %7.2f %7.2f   %5.2f%s\n",
           aoa, eL, eD, eD?eL/eD:0, pL, pD, pD?pL/pD:0, Cl,
           (aoa>=12 && aoa<=16)?"  <-stall":"");
  }
  mj_deleteModel(mw);

  // ---------- 2) Free glide: ellipsoid glider vs polar glider ----------
  printf("\nFree glide from 3 m, launched at 14 m/s:\n");

  mjModel* ge = mj_loadXML("/tmp/glider_ellipsoid.xml", NULL, err, sizeof(err));
  if(!ge){printf("load ellipsoid glider: %s\n",err);return 1;}
  double junk;
  double de = glide(ge, 0, &junk);
  printf("  ellipsoid model : glide distance %.1f m\n", de);
  mj_deleteModel(ge);

  mjModel* gp = mj_loadXML("/tmp/glider_polar.xml", NULL, err, sizeof(err));
  if(!gp){printf("load polar glider: %s\n",err);return 1;}
  // resolve aero geoms by name and precompute area/AR
  const char* names[]={"wing","tail","fin"};
  for (int i=0;i<3;i++){
    int g = mj_name2id(gp, mjOBJ_GEOM, names[i]);
    if(g<0) continue;
    aero_geom[n_aero]=g;
    aero_S[n_aero]=M_PI*gp->geom_size[3*g+0]*gp->geom_size[3*g+1];
    aero_AR[n_aero]=(2*gp->geom_size[3*g+1])/(2*gp->geom_size[3*g+0]);
    n_aero++;
  }
  double pk;
  double dp = glide(gp, 1, &pk);
  printf("  airfoil-polar   : glide distance %.1f m  (peak wing AoA %.1f deg)\n", dp, pk);
  mj_deleteModel(gp);
  return 0;
}
