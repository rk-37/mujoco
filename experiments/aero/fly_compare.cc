// Live A/B viewer: red glider uses MuJoCo's built-in ellipsoid fluid model,
// blue glider gets a real airfoil polar (with stall) via mjcb_passive.
#include <cstdio>
#include <cstring>
#include <cmath>
#include <GLFW/glfw3.h>
#include <mujoco/mujoco.h>

#define DEG (M_PI/180.0)

mjModel* m = NULL;
mjData*  d = NULL;
mjvCamera cam; mjvOption opt; mjvScene scn; mjrContext con;
bool bl=false, bm=false, br=false; double lastx=0, lasty=0;
bool paused = false;

// --- airfoil polar (linear lift -> stall -> flat plate, + induced drag) ---
static void airfoil_polar(double a, double AR, double* Cl, double* Cd){
  const double as=14*DEG, Cla=2*M_PI, Cd0=0.02, e=0.9;
  double Ll=Cla*a, Dl=Cd0+Ll*Ll/(M_PI*AR*e);
  double Ls=2*sin(a)*cos(a), Ds=Cd0+2*sin(a)*sin(a);
  double s=1.0/(1.0+exp(-(fabs(a)-as)/(3*DEG)));
  *Cl=(1-s)*Ll+s*Ls; *Cd=(1-s)*Dl+s*Ds;
}
static int aero[4], nA=0; double aS[4], aAR[4];
static double lastAoA=0, lastL=0;
static void polar_passive(const mjModel* m, mjData* d){
  double rho=m->opt.density;
  for(int i=0;i<nA;i++){
    int g=aero[i]; const mjtNum* R=d->geom_xmat+9*g;
    double cx[3]={R[0],R[3],R[6]}, sy[3]={R[1],R[4],R[7]}, nz[3]={R[2],R[5],R[8]};
    mjtNum v6[6]; mj_objectVelocity(m,d,mjOBJ_GEOM,g,v6,0);
    double Vr[3]={v6[3]-m->opt.wind[0],v6[4]-m->opt.wind[1],v6[5]-m->opt.wind[2]};
    double V=sqrt(Vr[0]*Vr[0]+Vr[1]*Vr[1]+Vr[2]*Vr[2]); if(V<1e-4) continue;
    double vc=Vr[0]*cx[0]+Vr[1]*cx[1]+Vr[2]*cx[2];
    double vn=Vr[0]*nz[0]+Vr[1]*nz[1]+Vr[2]*nz[2];
    double a=atan2(-vn,vc), Cl,Cd; airfoil_polar(a,aAR[i],&Cl,&Cd);
    double q=0.5*rho*V*V, S=aS[i], L=q*S*Cl, D=q*S*Cd;
    double dd[3]={-Vr[0]/V,-Vr[1]/V,-Vr[2]/V}, ld[3]; mju_cross(ld,Vr,sy);
    double ln=sqrt(ld[0]*ld[0]+ld[1]*ld[1]+ld[2]*ld[2]); if(ln<1e-6) continue;
    for(int k=0;k<3;k++) ld[k]/=ln;
    mjtNum f[3]={dd[0]*D+ld[0]*L,dd[1]*D+ld[1]*L,dd[2]*D+ld[2]*L}, t[3]={0,0,0};
    mj_applyFT(m,d,f,t,d->geom_xpos+3*g,m->geom_bodyid[g],d->qfrc_passive);
    if(i==0){ lastAoA=a/DEG; lastL=L; }   // remember main-wing state for HUD
  }
}

void keyboard(GLFWwindow* w,int key,int sc,int act,int mods){
  if(act!=GLFW_PRESS) return;
  if(key==GLFW_KEY_BACKSPACE){ mj_resetDataKeyframe(m,d,0); mj_forward(m,d); }
  if(key==GLFW_KEY_SPACE) paused=!paused;
}
void mouse_button(GLFWwindow* w,int b,int act,int mods){
  bl=glfwGetMouseButton(w,GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS;
  bm=glfwGetMouseButton(w,GLFW_MOUSE_BUTTON_MIDDLE)==GLFW_PRESS;
  br=glfwGetMouseButton(w,GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS;
  glfwGetCursorPos(w,&lastx,&lasty);
}
void mouse_move(GLFWwindow* w,double x,double y){
  if(!bl&&!bm&&!br) return;
  double dx=x-lastx, dy=y-lasty; lastx=x; lasty=y;
  int W,H; glfwGetWindowSize(w,&W,&H);
  bool sh=glfwGetKey(w,GLFW_KEY_LEFT_SHIFT)==GLFW_PRESS||glfwGetKey(w,GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS;
  mjtMouse a = br?(sh?mjMOUSE_MOVE_H:mjMOUSE_MOVE_V):bl?(sh?mjMOUSE_ROTATE_H:mjMOUSE_ROTATE_V):mjMOUSE_ZOOM;
  mjv_moveCamera(m,a,dx/H,dy/H,&scn,&cam);
}
void scroll(GLFWwindow* w,double xo,double yo){ mjv_moveCamera(m,mjMOUSE_ZOOM,0,-0.05*yo,&scn,&cam); }

int main(){
  char err[1000]="";
  m=mj_loadXML("/tmp/fly_compare.xml",0,err,1000);
  if(!m) mju_error("load: %s",err);
  d=mj_makeData(m);

  // resolve polar aero geoms + body ids
  const char* nm[]={"wing","tail","fin"};
  for(int i=0;i<3;i++){ int g=mj_name2id(m,mjOBJ_GEOM,nm[i]); if(g<0)continue;
    aero[nA]=g; aS[nA]=M_PI*m->geom_size[3*g]*m->geom_size[3*g+1];
    aAR[nA]=m->geom_size[3*g+1]/m->geom_size[3*g]; nA++; }
  int be=mj_name2id(m,mjOBJ_BODY,"glider_e");
  int bp=mj_name2id(m,mjOBJ_BODY,"glider_p");
  mjcb_passive = polar_passive;

  mj_resetDataKeyframe(m,d,0);

  if(!glfwInit()) mju_error("glfw init");
  GLFWwindow* win=glfwCreateWindow(1280,800,
    "Ellipsoid (red) vs Airfoil polar (blue)  -  SPACE pause, BACKSPACE relaunch",NULL,NULL);
  glfwMakeContextCurrent(win); glfwSwapInterval(1);
  mjv_defaultCamera(&cam); mjv_defaultOption(&opt);
  mjv_defaultScene(&scn); mjr_defaultContext(&con);
  mjv_makeScene(m,&scn,3000); mjr_makeContext(m,&con,mjFONTSCALE_150);
  // chase camera following the blue (polar) glider
  cam.type=mjCAMERA_TRACKING; cam.trackbodyid=bp;
  cam.distance=6.0; cam.azimuth=90; cam.elevation=-12;
  glfwSetKeyCallback(win,keyboard); glfwSetCursorPosCallback(win,mouse_move);
  glfwSetMouseButtonCallback(win,mouse_button); glfwSetScrollCallback(win,scroll);

  double landT=-1;
  while(!glfwWindowShouldClose(win)){
    mjtNum t0=d->time;
    if(!paused){
      while(d->time-t0<1.0/60.0) mj_step(m,d);
      // auto-relaunch ~2 s after the blue glider lands
      if(d->xpos[3*bp+2]<0.12){ if(landT<0) landT=d->time;
        else if(d->time-landT>2.0){ mj_resetDataKeyframe(m,d,0); landT=-1; } }
      else landT=-1;
    }

    mjrRect vp={0,0,0,0}; glfwGetFramebufferSize(win,&vp.width,&vp.height);
    mjv_updateScene(m,d,&opt,NULL,&cam,mjCAT_ALL,&scn);
    mjr_render(vp,&scn,&con);

    // HUD
    double xe=d->xpos[3*be], ze=d->xpos[3*be+2];
    double xp=d->xpos[3*bp], zp=d->xpos[3*bp+2];
    int ae=m->body_dofadr[be], ap=m->body_dofadr[bp];
    double se=sqrt(d->qvel[ae]*d->qvel[ae]+d->qvel[ae+2]*d->qvel[ae+2]);
    double sp=sqrt(d->qvel[ap]*d->qvel[ap]+d->qvel[ap+2]*d->qvel[ap+2]);
    char l[512], r[256];
    std::snprintf(l,sizeof l,
      "RED  ellipsoid model\n  dist %.1f m   alt %.2f m   speed %.1f m/s\n"
      "BLUE airfoil polar (w/ stall)\n  dist %.1f m   alt %.2f m   speed %.1f m/s\n  wing AoA %.1f deg",
      xe,ze,se, xp,zp,sp, lastAoA);
    std::snprintf(r,sizeof r,"%s", paused?"PAUSED":"running");
    mjr_overlay(mjFONT_NORMAL,mjGRID_TOPLEFT,vp,l,0,&con);
    mjr_overlay(mjFONT_NORMAL,mjGRID_TOPRIGHT,vp,r,0,&con);

    glfwSwapBuffers(win); glfwPollEvents();
  }
  mjv_freeScene(&scn); mjr_freeContext(&con);
  mj_deleteData(d); mj_deleteModel(m);
#if defined(__APPLE__)||defined(_WIN32)
  glfwTerminate();
#endif
  return 0;
}
