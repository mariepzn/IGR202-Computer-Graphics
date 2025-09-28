#ifndef _FLUID_HPP_
#define _FLUID_HPP_

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <set>

#ifndef M_PI
    #define M_PI 3.141592
#endif

#include "Vector.hpp"
#include "RigidSolver.hpp"
class SphSolver {
public:
  explicit SphSolver(
    const Real nu=0.08, const Real h=0.5, const Real density=1e3,
    const Vec2f g=Vec2f(0, -9.8), const Real eta=0.01, const Real gamma=7.0) :
    _kernel(h), _nu(nu), _h(h), _d0(density),
    _g(g), _eta(eta), _gamma(gamma)
  {
    _dt = 0.001;
    _m0 = _d0*_h*_h;
    _c = std::fabs(_g.y)/_eta;
    _k = _d0*_c*_c/_gamma;
  }

  // assume an arbitrary grid with the size of res_x*res_y; a fluid mass fill up
  // the size of f_width, f_height; each cell is sampled with 2x2 particles.
  void initScene(
    const int res_x, const int res_y, const int f_width, const int f_height, RigidSolver *object)
  {
    _pos.clear();

    _resX = res_x;
    _resY = res_y;
    
    
    // set wall for boundary
    _l = 0.5*_h;
    _r = static_cast<Real>(res_x) - 0.5*_h;
    _b = 0.5*_h;
    _t = static_cast<Real>(res_y) - 0.5*_h;

    // sample a fluid mass
    for(int j=0; j<f_height+4; ++j) {
      for(int i=0; i<res_x; ++i) {
        _pos.push_back(Vec2f(i+0.25, j+0.25));
        _pos.push_back(Vec2f(i+0.75, j+0.25));
        _pos.push_back(Vec2f(i+0.25, j+0.75));
        _pos.push_back(Vec2f(i+0.75, j+0.75));
      }
    }

    // make sure for the other particle quantities
    _vel = std::vector<Vec2f>(_pos.size(), Vec2f(0, 0));
    _acc = std::vector<Vec2f>(_pos.size(), Vec2f(0, 0));
   
    _p   = std::vector<Real>(_pos.size(), 0);
    _d   = std::vector<Real>(_pos.size(), 0);

    _col = std::vector<float>(_pos.size()*4, 1.0); // RGBA
    _vln = std::vector<float>(_pos.size()*4, 0.0); // GL_LINES
    

    updateColor();
    _object = object;
    _p0 = _object->body->pos2d;
    ObjectContrib=std::vector<float>(_object->getPos().size(), 0.0);

    
    ObjectContrib=_object->getContribution();
    
    //_forceObject = std::vector<Vec2f>(_p0.size(), Vec2f(0.0));
    buildNeighborCellsFluid();
  }
   void update(bool gShowVel)
  { 

  _forceObject = std::vector<Vec2f>(_p0.size(), Vec2f(0.0));
  _object->buildCellsSolid();

  _acc = std::vector<Vec2f>(_pos.size(), Vec2f(0, 0));
  
   _neighbors_Fluid = std::vector<std::vector<tIndex>>(_pos.size(), std::vector<tIndex>(0));
   _neighbors_Object = std::vector<std::vector<tIndex>>(_pos.size(), std::vector<tIndex>(0));
 
  std::map<tIndex, std::set<tIndex>>  object = _object->getpidXinGridSolid();
  std::vector<Vec2f>  posBody = _object->getPos();


   for(int i=0 ; i<_pos.size();i++){
    int x = floor(_pos[i].x);
    int y = floor(_pos[i].y);
    std::vector<tIndex> cellsNeighbors = getNeighbors(x,y);
      for(auto &n:cellsNeighbors){
        std::set<tIndex> object = _object->getpidXinGridSolid()[n];
        std::vector<Vec2f>  posBody = _object->getPos();
        std::cout<<posBody[1]<<std::endl;
        for(auto &p :_pidxInGridFluid[n]){
          if((_pos[i]-_pos[p]).length()<=_kernel.supportRadius()){
        _neighbors_Fluid[i].push_back(p);}}
        for(auto p :object ){
          if((_pos[i]-posBody[p]).length()<=_kernel.supportRadius()){
          _neighbors_Object[i].push_back(p);
}}}} 
  //std::cout<<_neighbors_Fluid.size()<<std::endl;

    computeDensity();
    
    computePressure();  
    applyBodyForce();
    applyPressureForce();
    applyViscousForce();
    updateVelocity();
    updateCell();
    resolveCollision();

    updateColor();
    if(gShowVel) updateVelLine();
  }

tIndex particleCount() const { return _pos.size(); }
const Vec2f& position(const tIndex i) const { return _pos[i]; }
const float& color(const tIndex i) const { return _col[i]; }
const float& vline(const tIndex i) const { return _vln[i]; }

int resX() const { return _resX; }
int resY() const { return _resY; }

    void buildNeighborCellsFluid()
  {
    {
    //_pidxInGridFluid.clear();
    //assign the cell x to the particle i
    for(tIndex i=0; i<_pos.size();++i){
      tIndex x = idx1d(floor(_pos[i].x), floor(_pos[i].y));
      _pidxInGridFluid[x].insert(i);
    }
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
      if(y<_resY-1){
        neighbors.push_back(idx1d(x-1,y+1)); //top left corner
      }
    }
    if (x<_resX-1){
      neighbors.push_back(idx1d(x+1,y)); //right cell
      if (y>0){
        neighbors.push_back(idx1d(x+1,y-1)); //lower right corner
      }
      if(y<(_resY-1)){
        neighbors.push_back(idx1d(x+1,y+1)); //top right corner
      }
    }
    if(y>0){
      neighbors.push_back(idx1d(x,y-1));//bottom cell
    }
    if(y<(_resY-1)){
      neighbors.push_back(idx1d(x,y+1));//top cell
    }
    return neighbors;
  }

    void computeDensity()
  {
    
    for (tIndex i =0; i<_pos.size(); i++){
      _d[i] = 0.0;
        for(int p : _neighbors_Fluid[i]){
          _d[i] += _m0* _kernel.w(_pos[i]-_pos[p]);
        }
        for(int p : _neighbors_Object[i]){
          _d[i] += ObjectContrib[p] * _kernel.w(_pos[i]-_object->getPos()[p]);
          
        }
      
    }
  }

  void computePressure()
  {
    for(tIndex p = 0; p < _pos.size(); ++p){
      _p[p]= std::max(0.0f, _k*(pow(_d[p]/_d0, _gamma)-1.0f));
    }
  }

  void applyBodyForce()
  {
    for(tIndex p = 0; p < _pos.size(); ++p){
      _acc[p] += _g;
    }
  }

  void applyPressureForce()
  {
    for (tIndex i = 0; i < _pos.size(); i++) {
        for (int p : _neighbors_Fluid[i]) {
          if (i != p)
            _acc[i] -= _m0 * ( _p[i] / (_d[i] * _d[i]) + _p[p] / (_d[p] * _d[p]))
                   * _kernel.grad_w(_pos[i] - _pos[p]);
        }
        for (int p : _neighbors_Object[i]){
          
            Vec2f F = _object->getContribution()[p]*_p[i]/(_d[i]*_d[i])*_kernel.grad_w(_pos[i]-_object->getPos()[p]);
            _acc[i]-=F;
            //std::cout<<"force"<<F[1]<<std::endl;
            _forceObject[p]+=_m0*F;
          }

    }


  }

  void applyViscousForce()
  {
   for (tIndex i = 0; i < _pos.size(); i++) {
      
      
        for (auto &p : _neighbors_Fluid[i]) {                           
          if (i != p)
          _acc[i] += 2.0 * _nu * _m0 / (_d[i] + _d[p])
                   * (_pos[i]-_pos[p]).dotProduct(_kernel.grad_w((_pos[i]-_pos[p])))
                   / ((_pos[i]-_pos[p]).dotProduct((_pos[i]-_pos[p])) + 0.01 * _h * _h)
                   * (_vel[i] - _vel[p]);
        }
        for(auto&p:_neighbors_Object[i]){

         Vec2f force = 2.0 * _nu / _d[i] * _object->getContribution()[p]
                   * (_pos[i]-_object->getPos()[p]).dotProduct(_kernel.grad_w((_pos[i]-_object->getPos()[p])))
                   / ((_pos[i]-_object->getPos()[p]).dotProduct((_pos[i]-_object->getPos()[p])) + 0.01 * _h * _h)
                   * (_vel[i] - _object->velocity(_object->getPos()[p]));

        _acc[i] += force;
        _forceObject[p]-=force*_m0;

        }
      
    }
  }

  void updateVelocity()
  {
   for (tIndex p = 0; p < _pos.size(); p++)
      _vel[p] = _vel[p] + _dt*_acc[p];
  }

  void updateCell()
  {
    std::vector<Vec2f> pos =_pos;

    for (tIndex p = 0; p < _pos.size(); p++){
      _pos[p] = _pos[p] + _dt*_vel[p];

    if(static_cast<int>(pos[p].x) != static_cast<int>(_pos[p].x)||static_cast<int>(pos[p].y) != static_cast<int>(_pos[p].y)){
      tIndex x = idx1d(floor(pos[p].x), floor(pos[p].y));
      tIndex newx = idx1d(floor(_pos[p].x), floor(_pos[p].y));
      _pidxInGridFluid[x].erase(p);
      _pidxInGridFluid[newx].insert(p);
    }
    }
  }
 std::vector<Vec2f> getForce(){
  return _forceObject;
 }

   // simple collision detection/resolution for each particle
  void resolveCollision()
  {
    std::vector<tIndex> need_res;
    for(tIndex i=0; i<particleCount(); ++i) {
      if(_pos[i].x<_l || _pos[i].y<_b || _pos[i].x>_r || _pos[i].y>_t)
        need_res.push_back(i);
    }

    for(
      std::vector<tIndex>::const_iterator it=need_res.begin();
      it<need_res.end();
      ++it) {
      const Vec2f p0 = _pos[*it];
      _pos[*it].x = clamp(_pos[*it].x, _l, _r);
      _pos[*it].y = clamp(_pos[*it].y, _b, _t);
      _vel[*it] = (_pos[*it] - p0)/_dt;
    }
  }

  void updateColor()
  {
    for(tIndex i=0; i<particleCount(); ++i) {
      _col[i*4+0] = 0.6;
      _col[i*4+1] = 0.6;
      _col[i*4+2] = _d[i]/_d0;
    }
  }
  Real get_dt() const { return _dt; }
  void updateVelLine()
  {
    for(tIndex i=0; i<particleCount(); ++i) {
      _vln[i*4+0] = _pos[i].x;
      _vln[i*4+1] = _pos[i].y;
      _vln[i*4+2] = _pos[i].x + _vel[i].x;
      _vln[i*4+3] = _pos[i].y + _vel[i].y;
    }
  }

  inline tIndex idx1d(const int i, const int j) { if (i < 0 || j < 0) {
      return -1;
    }
    if (i >= resX() || j >= resY()) {
      return -1;
    }

    return i + j * resX();}

const CubicSpline _kernel; //Kernel
RigidSolver *_object; //Rigid Object 
std::vector<Vec2f> _forceObject; //Forces applied to the object
std::vector<Vec2f> _p0; //Initial position of the solid
std::vector<std::vector<tIndex>>  _neighbors_Fluid; //Fluid Neighbors of the fluid particles
std::vector<std::vector<tIndex>>  _neighbors_Object; //Rigid neighbors of the fluid particles

  // particle data
  std::vector<Vec2f> _pos;     // position
  std::vector<Vec2f> _vel;      // velocity
  std::vector<Vec2f> _acc;      // acceleration
  std::vector<Real>  _p;        // pressure
  std::vector<Real>  _d;        // density

   std::map<tIndex, std::set<tIndex>> _pidxInGridFluid; // will help find neighbor particles

  std::vector<float> _col;    // particle color; just for visualization
  std::vector<float> _vln;    // particle velocity lines; just for visualization

  // simulation
  Real _dt;                     // time step

  int _resX, _resY;             // background grid resolution

  // wall
  Real _l, _r, _b, _t;          // wall (boundary)

  // SPH coefficients
  Real _nu;                     // viscosity coefficient
  Real _d0;                     // rest density
  Real _h;                      // particle spacing (i.e., diameter)
  Vec2f _g;                     // gravity

  Real _m0;                     // rest mass
  Real _k;                      // EOS coefficient

  Real _eta;
  Real _c;                      // speed of sound
  Real _gamma;                  // EOS power factor
  std::vector<float> ObjectContrib; //Volume contribution of the boundary particles 
};
#endif