// Target drone: a quadrotor cruising sideways at ~30 m, ~15 m/s, in steady state.
// Drag is present (inertia-box fluid model), so constant-velocity cruise needs a
// trimmed forward tilt: thrust T at pitch theta gives T*cos(theta)=mg (hold alt)
// and T*sin(theta)=drag(v) (hold speed). theta,T are found by a one-time trim solve
// at startup; flight then holds the 4 rotor controls CONSTANT (constant F).
//
// Run with no args for the GLFW viewer; run with "--headless" for a 10 s text check.
#include <cstdio>
#include <cstring>
#include <cmath>
#include <GLFW/glfw3.h>
#include <mujoco/mujoco.h>

#define DEG (M_PI/180.0)
static const char* MODEL = "/Users/rihards/src/mujoco/experiments/interceptor/target.xml";

mjModel* m = NULL;
mjData*  d = NULL;
mjvCamera cam; mjvOption opt; mjvScene scn; mjrContext con;
bool bl=false, bm=false, br=false; double lastx=0, lasty=0; bool paused=false;

int    g_bid=0, g_adr=0;          // target body id and its qvel/qacc base index
double g_theta=0, g_Tper=0;       // solved trim pitch (rad) and per-rotor thrust (N)
double g_v=15.0, g_alt=30.0;      // cruise speed and altitude
double g_x0=-150.0, g_xend=150.0; // start x (far approach) and reset x (far departure)
double g_yoff=30.0;               // lateral offset of the flight line (so not over zenith)
double g_obs[3]={0,0,1.7};        // ground observer eye position (on the world surface)

// instantaneous world-frame linear accel for a candidate pitch. For this theta we
// pick total thrust T so vertical accel is ~0 (accounts for the small vertical fluid
// force at nonzero pitch), then return the residual horizontal accel.
static void eval_accel(double theta, double mg, double* ax, double* az, double* Tper){
  double mass = mj_getTotalmass(m);
  double T = mg/cos(theta);                 // first guess: pure vertical balance
  d->qpos[0]=0; d->qpos[1]=0; d->qpos[2]=g_alt;
  d->qpos[3]=cos(theta/2); d->qpos[4]=0; d->qpos[5]=sin(theta/2); d->qpos[6]=0;
  d->qvel[0]=g_v; d->qvel[1]=0; d->qvel[2]=0;
  d->qvel[3]=0;   d->qvel[4]=0; d->qvel[5]=0;
  // Newton on T to zero vertical accel: d(az)/dT = cos(theta)/mass (linear => 1 step exact)
  for(int k=0;k<3;k++){
    for(int i=0;i<4;i++) d->ctrl[i]=T/4;
    mj_forward(m,d);
    double az_k = d->qacc[g_adr+2];
    if(fabs(az_k) < 1e-12) break;
    T -= az_k*mass/cos(theta);
  }
  *Tper = T/4;
  *ax = d->qacc[g_adr+0];
  *az = d->qacc[g_adr+2];
}

// bisection on horizontal accel -> (g_theta, g_Tper)
static void trim(){
  double mg = mj_getTotalmass(m) * (-m->opt.gravity[2]);
  // body-z = [sin(theta),0,cos(theta)] in world, so thrust pushes +x for theta>0.
  // Bracket straddles 0: ax(-25 deg)<0 (thrust -x + drag), ax(+25 deg)>0.
  double lo=-25*DEG, hi=25*DEG, ax, az, Tper, axlo, axhi;
  eval_accel(lo,mg,&axlo,&az,&Tper);
  eval_accel(hi,mg,&axhi,&az,&Tper);
  for(int it=0; it<60; it++){
    double mid=0.5*(lo+hi);
    eval_accel(mid,mg,&ax,&az,&Tper);
    if((ax>0)==(axlo>0)){ lo=mid; axlo=ax; } else { hi=mid; }
  }
  double mid=0.5*(lo+hi);
  eval_accel(mid,mg,&ax,&az,&Tper);
  g_theta=mid; g_Tper=Tper;
  printf("[trim] mass=%.3f kg  mg=%.3f N  theta=%.3f deg  T/rotor=%.4f N  resid ax=%.2e az=%.2e\n",
         mj_getTotalmass(m), mg, g_theta/DEG, g_Tper, ax, az);
}

// place the drone at the solved trim and load constant controls
static void apply_trim_state(){
  d->qpos[0]=g_x0; d->qpos[1]=g_yoff; d->qpos[2]=g_alt;
  d->qpos[3]=cos(g_theta/2); d->qpos[4]=0; d->qpos[5]=sin(g_theta/2); d->qpos[6]=0;
  d->qvel[0]=g_v; d->qvel[1]=0; d->qvel[2]=0;
  d->qvel[3]=0;   d->qvel[4]=0; d->qvel[5]=0;
  for(int i=0;i<4;i++) d->ctrl[i]=g_Tper;
  d->time=0;
  mj_forward(m,d);
}

// ---- headless 10 s check ----
static int run_headless(){
  apply_trim_state();
  // full traverse: -150 m -> +150 m at 15 m/s = 300 m / 20 s
  printf("\n   t       x      alt    speed   pitch(deg)   vz\n");
  double next=0, t=0;
  while(t <= 20.0001){
    if(t >= next-1e-9){
      double vx=d->qvel[0],vy=d->qvel[1],vz=d->qvel[2];
      double speed=sqrt(vx*vx+vy*vy+vz*vz);
      double pitch=2.0*atan2(d->qpos[5],d->qpos[3])/DEG;
      printf(" %5.1f  %7.1f  %6.2f  %6.2f   %+8.2f   %+6.3f\n",
             t, d->qpos[0], d->qpos[2], speed, pitch, vz);
      next += 1.0;
    }
    for(int i=0;i<4;i++) d->ctrl[i]=g_Tper;
    mj_step(m,d);
    t += m->opt.timestep;   // independent clock (d->time would reset on traverse loop)
  }
  return 0;
}

// ---- input callbacks ----
void keyboard(GLFWwindow* w,int key,int sc,int act,int mods){
  if(act!=GLFW_PRESS) return;
  if(key==GLFW_KEY_SPACE) paused=!paused;
  if(key==GLFW_KEY_BACKSPACE) apply_trim_state();
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

int main(int argc, char** argv){
  char err[1000]="";
  m=mj_loadXML(MODEL,0,err,1000);
  if(!m) mju_error("load: %s",err);
  d=mj_makeData(m);
  g_bid=mj_name2id(m,mjOBJ_BODY,"target");
  g_adr=m->body_dofadr[g_bid];

  trim();

  if(argc>1 && !strcmp(argv[1],"--headless")){
    int rc=run_headless();
    mj_deleteData(d); mj_deleteModel(m);
    return rc;
  }

  apply_trim_state();

  if(!glfwInit()) mju_error("glfw init");
  GLFWwindow* win=glfwCreateWindow(1280,800,
    "Target drone  -  SPACE pause, BACKSPACE relaunch at trim",NULL,NULL);
  glfwMakeContextCurrent(win); glfwSwapInterval(1);
  mjv_defaultCamera(&cam); mjv_defaultOption(&opt);
  mjv_defaultScene(&scn); mjr_defaultContext(&con);
  mjv_makeScene(m,&scn,3000); mjr_makeContext(m,&con,mjFONTSCALE_150);

  cam.type=mjCAMERA_FREE;   // re-aimed each frame from the ground observer (below)

  glfwSetKeyCallback(win,keyboard); glfwSetCursorPosCallback(win,mouse_move);
  glfwSetMouseButtonCallback(win,mouse_button); glfwSetScrollCallback(win,scroll);

  while(!glfwWindowShouldClose(win)){
    mjtNum t0=d->time;
    if(!paused){
      for(int i=0;i<4;i++) d->ctrl[i]=g_Tper;   // hold constant thrust
      while(d->time-t0<1.0/60.0) mj_step(m,d);
      if(d->qpos[0] >= g_xend) apply_trim_state();   // far departure -> reset & loop
    }

    // ground-observer camera: eye planted on the surface, gaze follows the drone
    double D[3]={d->xpos[3*g_bid+0],d->xpos[3*g_bid+1],d->xpos[3*g_bid+2]};
    double vv[3]={D[0]-g_obs[0],D[1]-g_obs[1],D[2]-g_obs[2]};
    double range=sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
    double elev=asin(vv[2]/range)/DEG;
    cam.lookat[0]=D[0]; cam.lookat[1]=D[1]; cam.lookat[2]=D[2];
    cam.distance=range;
    cam.azimuth=atan2(vv[1],vv[0])/DEG;
    cam.elevation=elev;

    mjrRect vp={0,0,0,0}; glfwGetFramebufferSize(win,&vp.width,&vp.height);
    mjv_updateScene(m,d,&opt,NULL,&cam,mjCAT_ALL,&scn);
    mjr_render(vp,&scn,&con);

    double x=d->qpos[0], y=d->qpos[1], z=d->qpos[2];
    double vx=d->qvel[0], vy=d->qvel[1], vz=d->qvel[2];
    double speed=sqrt(vx*vx+vy*vy+vz*vz);
    double pitch=2.0*atan2(d->qpos[5],d->qpos[3])/DEG;
    char l[512], rr[128];
    std::snprintf(l,sizeof l,
      "TARGET DRONE (constant-thrust trim)\n"
      "  pos    %7.1f %6.1f %6.1f m\n"
      "  alt    %6.2f m   speed %5.2f m/s\n"
      "  pitch  %+5.2f deg (trim %+5.2f)\n"
      "  thrust %.3f %.3f %.3f %.3f N\n"
      "  observer: range %5.0f m  elev %+4.0f deg",
      x,y,z, z,speed, pitch, g_theta/DEG,
      d->ctrl[0],d->ctrl[1],d->ctrl[2],d->ctrl[3],
      range, elev);
    std::snprintf(rr,sizeof rr,"%s", paused?"PAUSED":"running");
    mjr_overlay(mjFONT_NORMAL,mjGRID_TOPLEFT,vp,l,0,&con);
    mjr_overlay(mjFONT_NORMAL,mjGRID_TOPRIGHT,vp,rr,0,&con);

    glfwSwapBuffers(win); glfwPollEvents();
  }
  mjv_freeScene(&scn); mjr_freeContext(&con);
  mj_deleteData(d); mj_deleteModel(m);
#if defined(__APPLE__)||defined(_WIN32)
  glfwTerminate();
#endif
  return 0;
}
