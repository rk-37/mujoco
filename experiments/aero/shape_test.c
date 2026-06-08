// How does wing shape change behavior under each aero model?
#include <stdio.h>
#include <math.h>
#include <mujoco/mujoco.h>
#define DEG (M_PI/180.0)

static void airfoil_polar(double a,double AR,double*Cl,double*Cd){
  const double as=14*DEG,Cla=2*M_PI,Cd0=0.02,e=0.9;
  double Ll=Cla*a,Dl=Cd0+Ll*Ll/(M_PI*AR*e);
  double Ls=2*sin(a)*cos(a),Ds=Cd0+2*sin(a)*sin(a);
  double s=1.0/(1.0+exp(-(fabs(a)-as)/(3*DEG)));
  *Cl=(1-s)*Ll+s*Ls; *Cd=(1-s)*Dl+s*Ds;
}

// write a single-wing model with given semi-axes, measure ellipsoid forces
static void ellipsoid_force(double sx,double sy,double sz,double aoa,
                            double* L,double* D){
  char xml[1024];
  snprintf(xml,sizeof xml,
    "<mujoco><option density='1.225' viscosity='1.8e-5' gravity='0 0 0'"
    " integrator='implicitfast' wind='15 0 0'/>"
    "<worldbody><body pos='0 0 1'><freejoint/>"
    "<geom type='ellipsoid' size='%g %g %g' mass='0.2' fluidshape='ellipsoid'/>"
    "</body></worldbody></mujoco>", sx,sy,sz);
  FILE* f=fopen("/tmp/_w.xml","w"); fputs(xml,f); fclose(f);
  char err[500]=""; mjModel* m=mj_loadXML("/tmp/_w.xml",0,err,500);
  if(!m){printf("err %s\n",err);*L=*D=0;return;}
  mjData* d=mj_makeData(m);
  double a=aoa*DEG; d->qpos[3]=cos(a/2); d->qpos[5]=sin(a/2);
  mju_zero(d->qvel,m->nv); mj_forward(m,d);
  *D=d->qfrc_passive[0]; *L=d->qfrc_passive[2];
  mj_deleteData(d); mj_deleteModel(m);
}

int main(void){
  double aoa=6.0, q=0.5*1.225*15*15;   // 6 deg, 15 m/s
  struct { const char* name; double sx,sy,sz; } W[] = {
    {"baseline   AR 5.8", 0.12, 0.70, 0.008},
    {"high-AR    AR 23 ", 0.06, 1.40, 0.008},  // same area, long & thin
    {"low-AR     AR 1.5", 0.24, 0.35, 0.008},  // same area, stubby
    {"baseline THIN     ", 0.12, 0.70, 0.004},  // thin section
    {"baseline THICK    ", 0.12, 0.70, 0.050},  // fat section
  };
  printf("All wings ~same planform area (0.264 m2) except where noted. AoA=6 deg, V=15 m/s.\n\n");
  printf("                     |-- ellipsoid (built-in) --|--- airfoil polar ---|\n");
  printf("  wing                lift(N) drag(N)  L/D       lift(N) drag(N)  L/D\n");
  for(int i=0;i<5;i++){
    double eL,eD; ellipsoid_force(W[i].sx,W[i].sy,W[i].sz,aoa,&eL,&eD);
    double S=M_PI*W[i].sx*W[i].sy, AR=W[i].sy/W[i].sx;
    double Cl,Cd; airfoil_polar(aoa*DEG,AR,&Cl,&Cd);
    double pL=q*S*Cl, pD=q*S*Cd;
    printf("  %s  %7.2f %6.2f %6.2f     %7.2f %6.2f %6.2f\n",
           W[i].name, eL,eD, eD?eL/eD:0, pL,pD, pD?pL/pD:0);
  }
  printf("\nNote: 'THIN' vs 'THICK' change only the cross-section thickness (sz).\n");
  return 0;
}
