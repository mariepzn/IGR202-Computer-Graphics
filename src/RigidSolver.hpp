#ifndef _RIGIDSOLVER_HPP_
#define _RIGIDSOLVER_HPP_
#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include "Quaternion.hpp"
class CubicSpline
{
public:
  explicit CubicSpline(const Real h = 1) : _dim(2)
  {
    setSmoothingLen(h);
  }
  void setSmoothingLen(const Real h)
  {
    const Real h2 = square(h), h3 = h2 * h;
    _h = h;
    _sr = 2e0 * h;
    _c[0] = 2e0 / (3e0 * h);
    _c[1] = 10e0 / (7e0 * M_PI * h2);
    _c[2] = 1e0 / (M_PI * h3);
    _gc[0] = _c[0] / h;
    _gc[1] = _c[1] / h;
    _gc[2] = _c[2] / h;
  }
  Real smoothingLen() const { return _h; }
  Real supportRadius() const { return _sr; }

  Real f(const Real l) const
  {
    const Real q = l / _h;
    if (q < 1e0)
      return _c[_dim - 1] * (1e0 - 1.5 * square(q) + 0.75 * cube(q));
    else if (q < 2e0)
      return _c[_dim - 1] * (0.25 * cube(2e0 - q));
    return 0;
  }
  Real derivative_f(const Real l) const
  {
    const Real q = l / _h;
    if (q <= 1e0)
      return _gc[_dim - 1] * (-3e0 * q + 2.25 * square(q));
    else if (q < 2e0)
      return -_gc[_dim - 1] * 0.75 * square(2e0 - q);
    return 0;
  }

  Real w(const Vec2f &rij) const { return f(rij.length()); }
  Vec2f grad_w(const Vec2f &rij) const { return grad_w(rij, rij.length()); }
  Vec2f grad_w(const Vec2f &rij, const Real len) const
  {
    if (len < 1e-6) {
      return {0.0, 0.0};
    }

    return derivative_f(len) * rij / len;
  }

private:
  unsigned int _dim;
  Real _h, _sr, _c[3], _gc[3];
};
struct BodyAttributes {
  BodyAttributes() :
    X(0, 0, 0), R(Mat3f::I()), P(0, 0, 0), L(0, 0, 0),
    V(0, 0, 0), omega(0, 0, 0), F(0, 0, 0), tau(0, 0, 0), Q(1 , Vec3f(0, 0, 0)){}



  Real M;                      // mass
  Mat3f I0, I0inv;              // inertia tensor and its inverse in body space
  Mat3f Iinv;                   // inverse of inertia tensor

  // rigid body state
  Vec3f X;                      // position
  Mat3f R;                      // rotation
  Vec3f P;                      // linear momentum
  Vec3f L;                      // angular momentum
  Quaternion Q;
  // auxiliary quantities
  Vec3f V;                      // linear velocity
  Vec3f omega;                  // angular velocity

  // force and torque
  Vec3f F;                      // force
  Vec3f tau;                    // torque

  // mesh's vertices in body space
  
  std::vector<Vec2f> pos2d; //2D vector of the particles position
  std::vector<Vec3f> vdata0; //initial position
  std::vector<float> color;
  int _resX; //grid reoslution
  int _resY; 

};

class Box : public BodyAttributes {
public:
  explicit Box(
    const Real w=1.0, const Real h=1.0, const Real d=1.0, const Real mass =1000.0,
    const Vec3f v0=Vec3f(0, 0, 0), const Vec3f omega0=Vec3f(0, 0, 0)) :
    width(w), height(h), depth(d)
  {
      
    V = v0;                     // initial velocity
    omega = omega0;             // initial angular velocity
    M = mass;
    I0 =Mat3f(
      M*(2*M_PI*2*2)/12, 0, 0,
      0, M*(2*M_PI*2*2)/12, 0,
      0, 0, M*(2*M_PI*2*2)/12
    );
    I0inv =I0.inverse();
    Iinv = R * I0inv * R.transpose();

    const Real pi = M_PI;
    const Real angleStep = 2 * pi / 45; 

        for (int i = 0; i < 45; ++i) {
            Real angle = i * angleStep; 
            Real x =  2 * cos(angle);
            Real y =   2 * sin(angle);
            vdata0.push_back(Vec3f(x, y,0));
            
        }

  color = std::vector<float>(vdata0.size()*4,1);
  for (int i =0; i<vdata0.size();++i){
    color[i * 4 + 0] = 1.0;
    color[i * 4 + 1] = 0.6;
    color[i * 4 + 2] = 0.6;
  }
  }

  // rigid body property
  Real width, height, depth;
};

class RigidSolver {
public:
  explicit RigidSolver(
    BodyAttributes *body0=nullptr, const Vec3f g=Vec3f(0, 0, 0),const Real kernel_value = 0.5) :
    body(body0), _g(g), _step(0), _sim_t(0), _kernel(kernel_value) {}

  void init(BodyAttributes *body0,  Vec3f center,int resX,int resY)
  {
    body = body0;
    body->_resX=resX;
    body->_resY=resY;
    
    
    
    body->pos2d = std::vector<Vec2f>(body->vdata0.size());
    for(int i=0; i<body->vdata0.size();i++){
      body->pos2d[i]= (body->R * body->vdata0[i] + body->X).get2d();
      //std::cout<<"R"<<body->X<<std::endl;
      }
    body->X = center;
    body-> _resX=resX;
    body->_resY=resY;

    body->F = Vec3f(0.0);
    body->tau = Vec3f(0.0f);

    _step = 0;
    _sim_t = 0;

    pos = std::vector<Vec2f>(body->vdata0.size(), Vec2f(0.0));
  
    buildCellsSolid();
    computeContribution();


  }
  inline tIndex idx1d(const int i, const int j) {
    if (i < 0 || j < 0) {
      return -1;
    }
    if (i >= body->_resX || j >= body->_resY) {
      return -1;
    }

    return i + j * body->_resX;
  }

  Vec2f velocity(Vec2f p) const {
    Vec3f p3d = Vec3f(p.x, p.y, 0.0);
    Vec3f v3d = body->V + (body->omega).crossProduct(p3d - body->X);
    return v3d.get2d();
  }

  Vec2f translate(Vec2f p) {
    Vec3f p3d = Vec3f(p.x, p.y, 0.0);
    p3d = body->R * p3d + body->X;
    return p3d.get2d();
  }
  void buildCellsSolid()
    {
  
    _pidxInGridSolid.clear();
    //assign the cell x to the particle i
    for(tIndex i=0; i<body->vdata0.size();++i){

      pos[i] = translate(body->pos2d[i]);
      //<<i<<std::endl;
      tIndex x = idx1d(floor(pos[i].x), floor(pos[i].y));
      
      _pidxInGridSolid[x].insert(i);
    }
  }

std::vector<tIndex> getNeighbors(int x, int y ){
    std::vector<tIndex> neighbors;
    neighbors.push_back(idx1d(x,y));
    if(x>0){
      neighbors.push_back(idx1d(x-1,y)); //left cell
      if (y>0){
        neighbors.push_back(idx1d(x-1,y-1)); //lower left corner
      }
      if(y<body->_resY-1){
        neighbors.push_back(idx1d(x-1,y+1)); //top left corner
      }
    }
    if (x<body->_resX-1){
      neighbors.push_back(idx1d(x+1,y)); //right cell
      if (y>0){
        neighbors.push_back(idx1d(x+1,y-1)); //lower right corner
      }
      if(y<(body->_resY-1)){
        neighbors.push_back(idx1d(x+1,y+1)); //top right corner
      }
    }
    if(y>0){
      neighbors.push_back(idx1d(x,y-1));//bottom cell
    }
    if(y<(body->_resY-1)){
      neighbors.push_back(idx1d(x,y+1));//top cell
    }
    return neighbors;
  }
  tIndex particleCount()  { return pos.size(); }
void computeContribution(){
  contribution = std::vector<float>(body->vdata0.size(), 0.0);
  for (int i = 0; i < body->vdata0.size(); ++i) {
    int x = floor(pos[i].x);
    int y = floor(pos[i].y);
    Real volume = 0.0;
    std::vector<tIndex> cellsNeighbors = getNeighbors(x,y);
    for(auto &n:cellsNeighbors){
      
      std::set<tIndex>neighbor_particles_Object=_pidxInGridSolid[n];
      
      for(auto &p:neighbor_particles_Object){
        
          volume += _kernel.w(pos[i] - pos[p]);
      }
    }
    contribution[i] = _d0 / volume;
      //std::cout<<"contrib"<<contribution[34]<<std::endl;
    }

}
  std::vector<Real> getContribution(){
    return contribution;
  }
std::map<tIndex, std::set<tIndex>>  getpidXinGridSolid(){
    return _pidxInGridSolid;
  }
  std::vector<Vec2f> getPos(){ return pos; }
  tIndex particleCount() const { return (tIndex)pos.size(); }
  void step(const Real dt,std::vector<Vec2f>  forces)
  {
    

    computeForceAndTorque(forces);



    body->P = body->P+dt*body->F;
    body->V = body->P / body->M;
    body->X = body->X + dt * body->V;
    

    body->Iinv = body->R * body->I0inv * body->R.transpose();
    body->L = body->L + dt * body->tau;
    body->omega = body->Iinv * body->L;
    body->Q += Quaternion(0, body->omega) * body->Q * dt / 2;
    body->R = body->Q.normalize().rotation();

    
    
    //body->R += body->omega.crossProductMatrix() * body->R * dt;

    ++_step;
    _sim_t += dt;

  }

BodyAttributes *body;

private:
  void computeForceAndTorque(std::vector<Vec2f> forces)
  { 

    body->F = Vec3f(0, 0, 0);
    body->F = body->M * _g;
    body->tau = Vec3f(0, 0, 0);
    for (int i = 0;i<body->vdata0.size();i++){
      Vec3f F = Vec3f(forces[i].x, forces[i].y, 0.0);
      Vec3f p = Vec3f(pos[i].x, pos[i].y, 0.0);
      body->F+=F;
      body->tau +=(p-body->X).crossProduct(F);
    }


    // instance force at the very first step
    /*if(_step==1) {
      body->F = Vec3f(0.15f, 0.15f, 0.03);
      body->tau = crossProductMatrix(body->R * body->vdata0[0]) * body->F;
    }*/
  }

  // simulation parameters
Vec3f _g;                     // gravity
tIndex _step;                 // step count
Real _sim_t;  
std::map<tIndex, std::set<tIndex>> _pidxInGridSolid; //Assing the particles ot each cell
const CubicSpline _kernel ; //kernel
std::vector<Vec2f> pos; //Actual position of the particles
std::vector<Real> contribution; //Volume contribution of the particles, used in the SphSolver
const Real _d0 = 1000; //Fluid density 

};

#endif  /* _RIGIDSOLVER_HPP_ */
