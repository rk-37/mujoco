// Interceptor scenario.
//   - target drone cruises sideways at 30 m, 15 m/s (constant-thrust trim).
//   - a small throwable interceptor is launched from a ground observer at ~30 mph
//     in the (noisy) direction of the target. Its motors + onboard camera point
//     along +body-z, so the camera looks where it is thrown.
//   - PLACEHOLDER guidance: it flies in a fixed direction with constant thrust and
//     misses. homing_control() is left empty for the future FPV-based controller.
// View: ground-observer main view + a picture-in-picture of the drone's up-camera
//       (bottom-right). Scenery (pylons/clouds/hills) gives motion reference.
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <GLFW/glfw3.h>
#include <mujoco/mujoco.h>

#define DEG (M_PI/180.0)
static const char* MODEL = "/Users/rihards/src/mujoco/experiments/interceptor/intercept.xml";

mjModel* m=NULL; mjData* d=NULL;
mjvCamera cam; mjvOption opt; mjvScene scn, scnf; mjrContext con;
bool bl=false, bm=false, br=false; double lastx=0, lasty=0; bool paused=false;

// target
int t_bid, t_qadr, t_vadr;
double g_theta=0, g_Tper=0, g_v=15.0, g_alt=30.0, g_x0=-150.0, g_xend=150.0, g_yoff=30.0;
// interceptor
int i_bid, i_qadr, i_vadr, i_act[4], i_camid;
double g_obs[3]={0,0,1.7};            // ground observer eye
double g_launch_from[3]={2,0,1.2};    // throw origin (beside the observer)
double g_throw_speed=13.41;           // 30 mph
double g_i_mg=0;                       // interceptor weight (set after load)
bool   i_launched=false;
double g_aim[3]={1,0,0};               // static (noisy) launch direction, fixed at throw

// ---------- small helpers ----------
static double frand_gauss(){           // ~N(0,1), Box-Muller
  double u1=(rand()+1.0)/(RAND_MAX+2.0), u2=(rand()+1.0)/(RAND_MAX+2.0);
  return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}
// quaternion (w,x,y,z) rotating +z onto unit vector dir
static void quat_z_to(const double dir[3], double q[4]){
  double dot=dir[2];
  if(dot> 0.999999){ q[0]=1;q[1]=q[2]=q[3]=0; return; }
  if(dot<-0.999999){ q[0]=0;q[1]=1;q[2]=q[3]=0; return; }   // 180 about x
  double ax[3]={-dir[1],dir[0],0};
  double an=sqrt(ax[0]*ax[0]+ax[1]*ax[1]); ax[0]/=an; ax[1]/=an;
  double ang=acos(dot), s=sin(ang/2);
  q[0]=cos(ang/2); q[1]=ax[0]*s; q[2]=ax[1]*s; q[3]=0;
}

// ---------- target trim (same scheme as target_sim.cc) ----------
static void eval_accel(double theta, double mg, double* ax, double* az, double* Tper){
  double mass = 1.2;                  // target mass (core+arms+rotors)
  double T = mg/cos(theta);
  d->qpos[t_qadr+0]=0; d->qpos[t_qadr+1]=0; d->qpos[t_qadr+2]=g_alt;
  d->qpos[t_qadr+3]=cos(theta/2); d->qpos[t_qadr+4]=0; d->qpos[t_qadr+5]=sin(theta/2); d->qpos[t_qadr+6]=0;
  d->qvel[t_vadr+0]=g_v; d->qvel[t_vadr+1]=0; d->qvel[t_vadr+2]=0;
  d->qvel[t_vadr+3]=0;   d->qvel[t_vadr+4]=0; d->qvel[t_vadr+5]=0;
  for(int k=0;k<3;k++){
    for(int i=0;i<4;i++) d->ctrl[i]=T/4;
    mj_forward(m,d);
    double az_k=d->qacc[t_vadr+2];
    if(fabs(az_k)<1e-12) break;
    T -= az_k*mass/cos(theta);
  }
  *Tper=T/4; *ax=d->qacc[t_vadr+0]; *az=d->qacc[t_vadr+2];
}
static void trim(){
  double mg = 1.2 * (-m->opt.gravity[2]);
  double lo=-25*DEG, hi=25*DEG, ax, az, Tper, axlo;
  eval_accel(lo,mg,&axlo,&az,&Tper);
  for(int it=0; it<60; it++){
    double mid=0.5*(lo+hi); eval_accel(mid,mg,&ax,&az,&Tper);
    if((ax>0)==(axlo>0)){ lo=mid; axlo=ax; } else { hi=mid; }
  }
  double mid=0.5*(lo+hi); eval_accel(mid,mg,&ax,&az,&Tper);
  g_theta=mid; g_Tper=Tper;
  printf("[trim] target theta=%.3f deg  T/rotor=%.4f N\n", g_theta/DEG, g_Tper);
}

// ---------- state setters ----------
static void set_target_start(){
  d->qpos[t_qadr+0]=g_x0; d->qpos[t_qadr+1]=g_yoff; d->qpos[t_qadr+2]=g_alt;
  d->qpos[t_qadr+3]=cos(g_theta/2); d->qpos[t_qadr+4]=0; d->qpos[t_qadr+5]=sin(g_theta/2); d->qpos[t_qadr+6]=0;
  d->qvel[t_vadr+0]=g_v; d->qvel[t_vadr+1]=0; d->qvel[t_vadr+2]=0;
  d->qvel[t_vadr+3]=0;   d->qvel[t_vadr+4]=0; d->qvel[t_vadr+5]=0;
}
static void stow_interceptor(){
  i_launched=false;
  d->qpos[i_qadr+0]=g_launch_from[0]; d->qpos[i_qadr+1]=g_launch_from[1]; d->qpos[i_qadr+2]=g_launch_from[2];
  d->qpos[i_qadr+3]=1; d->qpos[i_qadr+4]=0; d->qpos[i_qadr+5]=0; d->qpos[i_qadr+6]=0;
  for(int k=0;k<6;k++) d->qvel[i_vadr+k]=0;
  for(int i=0;i<4;i++) d->ctrl[i_act[i]]=0;
}
static void launch_interceptor(){
  // aim at the target's current position, with directional noise (~4 deg)
  double T[3]={d->xpos[3*t_bid+0],d->xpos[3*t_bid+1],d->xpos[3*t_bid+2]};
  double dir[3]={T[0]-g_launch_from[0],T[1]-g_launch_from[1],T[2]-g_launch_from[2]};
  double n=sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
  for(int k=0;k<3;k++) dir[k]/=n;
  double sigma=4*DEG;
  for(int k=0;k<3;k++) dir[k]+=sigma*frand_gauss();
  n=sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
  for(int k=0;k<3;k++) dir[k]/=n;
  for(int k=0;k<3;k++) g_aim[k]=dir[k];
  // keep the drone LEVEL: motors (and the up-camera) point straight up, so the
  // constant hover thrust exactly cancels gravity and the throw velocity carries
  // it in a true straight line (no ground crash). A future homing controller will
  // tilt the body to steer; this static placeholder does not.
  d->qpos[i_qadr+3]=1; d->qpos[i_qadr+4]=0; d->qpos[i_qadr+5]=0; d->qpos[i_qadr+6]=0;
  // throw velocity (toward the noisy target direction)
  d->qvel[i_vadr+0]=g_throw_speed*dir[0];
  d->qvel[i_vadr+1]=g_throw_speed*dir[1];
  d->qvel[i_vadr+2]=g_throw_speed*dir[2];
  d->qvel[i_vadr+3]=0; d->qvel[i_vadr+4]=0; d->qvel[i_vadr+5]=0;
  i_launched=true;
  printf("[launch] dir=(%.2f %.2f %.2f) toward target at (%.0f %.0f %.0f)\n",
         dir[0],dir[1],dir[2], T[0],T[1],T[2]);
}

// PLACEHOLDER guidance: fixed direction, constant thrust. No homing yet.
// TODO: replace with a controller that consumes the FPV camera stream and steers
//       the interceptor (re-orient body-z + modulate per-rotor thrust) to the target.
static void homing_control(){
  // constant hover-magnitude thrust along the (fixed) body-z = launch direction.
  double per = g_i_mg/4.0;
  for(int i=0;i<4;i++) d->ctrl[i_act[i]]=per;
}

static void scene_reset(){
  set_target_start();
  stow_interceptor();
  d->time=0;
  mj_forward(m,d);
}

// ---------- input ----------
void keyboard(GLFWwindow* w,int key,int sc,int act,int mods){
  if(act!=GLFW_PRESS) return;
  if(key==GLFW_KEY_SPACE) paused=!paused;
  if(key==GLFW_KEY_BACKSPACE) scene_reset();
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

// step the scenario forward one control tick (target + interceptor logic + reset)
static double g_min_miss=1e9;
static void control_and_step(){
  for(int i=0;i<4;i++) d->ctrl[i]=g_Tper;
  if(!i_launched){
    stow_interceptor();
    double tp[3]={d->xpos[3*t_bid+0],d->xpos[3*t_bid+1],d->xpos[3*t_bid+2]};
    double dx=tp[0]-g_launch_from[0],dy=tp[1]-g_launch_from[1],dz=tp[2]-g_launch_from[2];
    if(sqrt(dx*dx+dy*dy+dz*dz)<80.0) launch_interceptor();
  } else {
    homing_control();
    double mx=d->xpos[3*i_bid]-d->xpos[3*t_bid];
    double my=d->xpos[3*i_bid+1]-d->xpos[3*t_bid+1];
    double mz=d->xpos[3*i_bid+2]-d->xpos[3*t_bid+2];
    double md=sqrt(mx*mx+my*my+mz*mz);
    if(md<g_min_miss) g_min_miss=md;
  }
  mjtNum t0=d->time;
  while(d->time-t0<1.0/60.0) mj_step(m,d);
}

static int run_headless(){
  printf("\n   t     tgt_x   int_pos(x,y,z)        int_spd   miss\n");
  double next=0, tclock=0;
  int loops=0;
  while(loops<1){
    if(tclock>=next-1e-9){
      printf(" %5.1f  %6.1f   %6.1f %6.1f %6.1f   %6.2f   %6.1f\n",
        tclock, d->xpos[3*t_bid],
        d->xpos[3*i_bid],d->xpos[3*i_bid+1],d->xpos[3*i_bid+2],
        sqrt(d->qvel[i_vadr]*d->qvel[i_vadr]+d->qvel[i_vadr+1]*d->qvel[i_vadr+1]+d->qvel[i_vadr+2]*d->qvel[i_vadr+2]),
        i_launched ? sqrt(
          (d->xpos[3*i_bid]-d->xpos[3*t_bid])*(d->xpos[3*i_bid]-d->xpos[3*t_bid])+
          (d->xpos[3*i_bid+1]-d->xpos[3*t_bid+1])*(d->xpos[3*i_bid+1]-d->xpos[3*t_bid+1])+
          (d->xpos[3*i_bid+2]-d->xpos[3*t_bid+2])*(d->xpos[3*i_bid+2]-d->xpos[3*t_bid+2])) : -1.0);
      next+=0.5;
    }
    control_and_step();
    tclock+=1.0/60.0;
    if(d->qpos[t_qadr+0]>=g_xend){ scene_reset(); loops++; }
    if(tclock>30) break;
  }
  printf("\nclosest approach (miss distance): %.2f m\n", g_min_miss);
  return 0;
}

int main(int argc, char** argv){
  srand((unsigned)time(NULL));
  setbuf(stdout,NULL);
  char err[1000]="";
  m=mj_loadXML(MODEL,0,err,1000);
  if(!m) mju_error("load: %s",err);
  d=mj_makeData(m);

  t_bid=mj_name2id(m,mjOBJ_BODY,"target");
  i_bid=mj_name2id(m,mjOBJ_BODY,"interceptor");
  t_qadr=m->jnt_qposadr[m->body_jntadr[t_bid]]; t_vadr=m->body_dofadr[t_bid];
  i_qadr=m->jnt_qposadr[m->body_jntadr[i_bid]]; i_vadr=m->body_dofadr[i_bid];
  const char* an[4]={"i_m_px","i_m_nx","i_m_py","i_m_ny"};
  for(int i=0;i<4;i++) i_act[i]=mj_name2id(m,mjOBJ_ACTUATOR,an[i]);
  i_camid=mj_name2id(m,mjOBJ_CAMERA,"fpv");

  // interceptor mass for hover thrust
  g_i_mg = (m->body_mass[i_bid]) * (-m->opt.gravity[2]);
  printf("[init] interceptor mass=%.3f kg  hover=%.3f N\n", m->body_mass[i_bid], g_i_mg);

  trim();
  scene_reset();

  if(argc>1 && !strcmp(argv[1],"--headless")){
    int rc=run_headless();
    mj_deleteData(d); mj_deleteModel(m);
    return rc;
  }

  if(!glfwInit()) mju_error("glfw init");
  GLFWwindow* win=glfwCreateWindow(1280,800,
    "Interceptor  -  SPACE pause, BACKSPACE reset",NULL,NULL);
  glfwMakeContextCurrent(win); glfwSwapInterval(1);
  mjv_defaultCamera(&cam); mjv_defaultOption(&opt);
  mjv_defaultScene(&scn); mjv_defaultScene(&scnf); mjr_defaultContext(&con);
  mjv_makeScene(m,&scn,4000); mjv_makeScene(m,&scnf,4000);
  mjr_makeContext(m,&con,mjFONTSCALE_150);
  cam.type=mjCAMERA_FREE;   // re-aimed each frame from the ground observer

  glfwSetKeyCallback(win,keyboard); glfwSetCursorPosCallback(win,mouse_move);
  glfwSetMouseButtonCallback(win,mouse_button); glfwSetScrollCallback(win,scroll);

  // fpv camera object (fixed onboard camera)
  mjvCamera fcam; mjv_defaultCamera(&fcam);
  fcam.type=mjCAMERA_FIXED; fcam.fixedcamid=i_camid;

  while(!glfwWindowShouldClose(win)){
    if(!paused){
      control_and_step();
      if(d->qpos[t_qadr+0] >= g_xend) scene_reset();   // target departed -> loop
    }

    // ----- main view: ground observer, gaze follows the target -----
    double Tp[3]={d->xpos[3*t_bid+0],d->xpos[3*t_bid+1],d->xpos[3*t_bid+2]};
    double vv[3]={Tp[0]-g_obs[0],Tp[1]-g_obs[1],Tp[2]-g_obs[2]};
    double range=sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
    cam.lookat[0]=Tp[0]; cam.lookat[1]=Tp[1]; cam.lookat[2]=Tp[2];
    cam.distance=range;
    cam.azimuth=atan2(vv[1],vv[0])/DEG;
    cam.elevation=asin(vv[2]/range)/DEG;

    mjrRect full={0,0,0,0}; glfwGetFramebufferSize(win,&full.width,&full.height);
    mjv_updateScene(m,d,&opt,NULL,&cam,mjCAT_ALL,&scn);
    mjr_render(full,&scn,&con);

    // ----- HUD -----
    double isp = sqrt(d->qvel[i_vadr]*d->qvel[i_vadr]+d->qvel[i_vadr+1]*d->qvel[i_vadr+1]+d->qvel[i_vadr+2]*d->qvel[i_vadr+2]);
    double miss[3]={d->xpos[3*i_bid]-Tp[0], d->xpos[3*i_bid+1]-Tp[1], d->xpos[3*i_bid+2]-Tp[2]};
    double missd=sqrt(miss[0]*miss[0]+miss[1]*miss[1]+miss[2]*miss[2]);
    char l[512];
    std::snprintf(l,sizeof l,
      "TARGET   alt %.1f m  speed %.1f m/s  (observer range %.0f m)\n"
      "INTERCEPTOR  %s   speed %.1f m/s\n"
      "  distance to target: %.1f m\n"
      "guidance: PLACEHOLDER (straight line, const thrust) -- homing TODO",
      Tp[2], g_v, range,
      i_launched?"LAUNCHED":"armed (waiting for target <80 m)", isp,
      missd);
    mjr_overlay(mjFONT_NORMAL,mjGRID_TOPLEFT,full,l,0,&con);

    // ----- picture-in-picture: drone up-camera (bottom-right) -----
    int pw=full.width/4, ph=full.height/4, mg=12;
    mjrRect pip={full.width-pw-mg, mg, pw, ph};
    mjrRect bord={pip.left-3, pip.bottom-3, pip.width+6, pip.height+6};
    mjr_rectangle(bord, 1,1,1,1);                       // white border
    mjv_updateScene(m,d,&opt,NULL,&fcam,mjCAT_ALL,&scnf);
    mjr_render(pip,&scnf,&con);
    mjr_overlay(mjFONT_NORMAL,mjGRID_BOTTOMLEFT,pip,"FPV  (drone up-camera)",0,&con);

    glfwSwapBuffers(win); glfwPollEvents();
  }
  mjv_freeScene(&scn); mjv_freeScene(&scnf); mjr_freeContext(&con);
  mj_deleteData(d); mj_deleteModel(m);
#if defined(__APPLE__)||defined(_WIN32)
  glfwTerminate();
#endif
  return 0;
}
