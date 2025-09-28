// ----------------------------------------------------------------------------
// main.cpp
//
//  Created on: Fri Jan 22 20:45:07 2021
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: SPH simulator (DO NOT DISTRIBUTE!)
//
// Copyright 2021-2024 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
//
// Modified by Marie Pizzini on 2024-01-19
// ----------------------------------------------------------------------------

#define _USE_MATH_DEFINES

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <memory>
#ifndef M_PI
#define M_PI 3.141592
#endif

// #include "Vector.hpp"
#include "SphSolver.hpp"

// window parameters
GLFWwindow *gWindow = nullptr;
int gWindowWidth = 1024;
int gWindowHeight = 768;

// timer
float gAppTimer = 0.0;
float gAppTimerLastClockTime;
bool gAppTimerStoppedP = true;

// global options
bool gPause = true;
bool gSaveFile = false;
bool gShowGrid = true;
bool gShowVel = false;
int gSavedCnt = 0;

const int kViewScale = 30;

RigidSolver bodySolver = RigidSolver(nullptr, Vec3f(0, -9.8, 0));
std::shared_ptr<BodyAttributes> rigidAtt = nullptr;
SphSolver sphSolver(0.08, 0.5, 1e3, Vec2f(0, -9.8), 0.01, 7.0);

void printHelp()
{
  std::cout << "> Help:" << std::endl
            << "    Keyboard commands:" << std::endl
            << "    * H: print this help" << std::endl
            << "    * P: toggle simulation" << std::endl
            << "    * G: toggle grid rendering" << std::endl
            << "    * V: toggle velocity rendering" << std::endl
            << "    * S: save current frame into a file" << std::endl
            << "    * Q: quit the program" << std::endl;
}

// Executed each time the window is resized. Adjust the aspect ratio and the rendering viewport to the current window.
void windowSizeCallback(GLFWwindow *window, int width, int height)
{
  gWindowWidth = width;
  gWindowHeight = height;
  glViewport(0, 0, static_cast<GLint>(gWindowWidth), static_cast<GLint>(gWindowHeight));
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, sphSolver.resX(), 0, sphSolver.resY(), 0, 1);
}

// Executed each time a key is entered.
void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
  if (action == GLFW_PRESS && key == GLFW_KEY_H)
  {
    printHelp();
  }
  else if (action == GLFW_PRESS && key == GLFW_KEY_S)
  {
    gSaveFile = !gSaveFile;
  }
  else if (action == GLFW_PRESS && key == GLFW_KEY_G)
  {
    gShowGrid = !gShowGrid;
  }
  else if (action == GLFW_PRESS && key == GLFW_KEY_V)
  {
    gShowVel = !gShowVel;
  }
  else if (action == GLFW_PRESS && key == GLFW_KEY_P)
  {
    gAppTimerStoppedP = !gAppTimerStoppedP;
    if (!gAppTimerStoppedP)
      gAppTimerLastClockTime = static_cast<float>(glfwGetTime());
  }
  else if (action == GLFW_PRESS && key == GLFW_KEY_Q)
  {
    glfwSetWindowShouldClose(window, true);
  }
}

void initGLFW()
{
  // Initialize GLFW, the library responsible for window management
  if (!glfwInit())
  {
    std::cerr << "ERROR: Failed to init GLFW" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Before creating the window, set some option flags
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  // glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // only if requesting 3.0 or above
  // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_ANY_PROFILE); // for OpenGL below 3.2
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // Create the window
  gWindowWidth = sphSolver.resX() * kViewScale;
  gWindowHeight = sphSolver.resY() * kViewScale;
  gWindow = glfwCreateWindow(
      sphSolver.resX() * kViewScale, sphSolver.resY() * kViewScale,
      "Basic SPH Simulator", nullptr, nullptr);
  if (!gWindow)
  {
    std::cerr << "ERROR: Failed to open window" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
  }

  // Load the OpenGL context in the GLFW window
  glfwMakeContextCurrent(gWindow);

  // not mandatory for all, but MacOS X
  glfwGetFramebufferSize(gWindow, &gWindowWidth, &gWindowHeight);

  // Connect the callbacks for interactive control
  glfwSetWindowSizeCallback(gWindow, windowSizeCallback);
  glfwSetKeyCallback(gWindow, keyCallback);

  std::cout << "Window created: " << gWindowWidth << ", " << gWindowHeight << std::endl;
}

void clear();

void exitOnCriticalError(const std::string &message)
{
  std::cerr << "> [Critical error]" << message << std::endl;
  std::cerr << "> [Clearing resources]" << std::endl;
  clear();
  std::cerr << "> [Exit]" << std::endl;
  std::exit(EXIT_FAILURE);
}

void initOpenGL()
{
  // Load extensions for modern OpenGL
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    exitOnCriticalError("[Failed to initialize OpenGL context]");

  glDisable(GL_CULL_FACE);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  glViewport(0, 0, static_cast<GLint>(gWindowWidth), static_cast<GLint>(gWindowHeight));
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, sphSolver.resX(), 0, sphSolver.resY(), 0, 1);
}

void init()
{
  rigidAtt = std::make_shared<Box>(4.f, 1.f, 0.f, 15000.);
  bodySolver.init(rigidAtt.get(), {15., 15, 0.0}, 30,30);
  
  sphSolver.initScene(30, 30, 12, 7, &bodySolver);
  std::cout<<"fluidinit"<<std::endl;

  initGLFW(); // Windowing system
  initOpenGL();
}

void clear()
{
  glfwDestroyWindow(gWindow);
  glfwTerminate();
}

// The main rendering call
void render()
{
  glClearColor(.4f, .4f, .4f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // grid guides
  if (gShowGrid)
  {
    glBegin(GL_LINES);
    for (int i = 1; i < sphSolver.resX(); ++i)
    {
      glColor3f(0.3, 0.3, 0.3);
      glVertex2f(static_cast<Real>(i), 0.0);
      glColor3f(0.3, 0.3, 0.3);
      glVertex2f(static_cast<Real>(i), static_cast<Real>(sphSolver.resY()));
    }
    for (int j = 1; j < sphSolver.resY(); ++j)
    {
      glColor3f(0.3, 0.3, 0.3);
      glVertex2f(0.0, static_cast<Real>(j));
      glColor3f(0.3, 0.3, 0.3);
      glVertex2f(static_cast<Real>(sphSolver.resX()), static_cast<Real>(j));
    }
    glEnd();
  }

  // render particles
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  glPointSize(0.20f * kViewScale);

  glColorPointer(4, GL_FLOAT, 0, &sphSolver.color(0));
  glVertexPointer(2, GL_FLOAT, 0, &sphSolver.position(0));
  glDrawArrays(GL_POINTS, 0, sphSolver.particleCount());

  glDisableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_VERTEX_ARRAY);
  glColor4f(1,0,1,1);
  glVertexPointer(2, GL_FLOAT, sizeof(Vec2f), &bodySolver.getPos()[0].x); 
  glDrawArrays(GL_POINTS, 0, bodySolver.getPos().size());

  glDisableClientState(GL_VERTEX_ARRAY);



  // velocity
  if (gShowVel)
  {
    glColor4f(0.0f, 0.0f, 0.5f, 0.2f);

    glEnableClientState(GL_VERTEX_ARRAY);

    glVertexPointer(2, GL_FLOAT, 0, &sphSolver.vline(0));
    glDrawArrays(GL_LINES, 0, sphSolver.particleCount() * 2);

    glDisableClientState(GL_VERTEX_ARRAY);
  }

  if (gSaveFile)
  {
    std::stringstream fpath;
    fpath << "s" << std::setw(4) << std::setfill('0') << gSavedCnt++ << ".tga";

    std::cout << "Saving file " << fpath.str() << " ... " << std::flush;
    const short int w = gWindowWidth;
    const short int h = gWindowHeight;
    std::vector<int> buf(w * h * 3, 0);
    glReadPixels(0, 0, w, h, GL_BGR, GL_UNSIGNED_BYTE, &(buf[0]));

    FILE *out = fopen(fpath.str().c_str(), "wb");
    short TGAhead[] = {0, 2, 0, 0, 0, 0, w, h, 24};
    fwrite(&TGAhead, sizeof(TGAhead), 1, out);
    fwrite(&(buf[0]), 3 * w * h, 1, out);
    fclose(out);
    gSaveFile = false;

    std::cout << "Done" << std::endl;
  }
}

// Update any accessible variable based on the current time
void update(const float currentTime)
{
  if (!gAppTimerStoppedP)
  {
    // NOTE: When you want to use application's dt ...
    const float dt = currentTime - gAppTimerLastClockTime;
    gAppTimerLastClockTime = currentTime;
    gAppTimer += dt;

    // solve 10 steps
    // for (int i = 0; i < 10; ++i)
      sphSolver.update(gShowVel);
      auto forces = sphSolver.getForce();
      bodySolver.step(sphSolver.get_dt(), forces);
      // bodySolver.step(dt);
      // rigidMat = rigidAtt->worldMat();
  }
}

int main(int argc, char **argv)
{
  init();
  while (!glfwWindowShouldClose(gWindow))
  {
    update(static_cast<float>(glfwGetTime()));
    render();
    glfwSwapBuffers(gWindow);
    glfwPollEvents();
  }
  clear();
  std::cout << " > Quit" << std::endl;
  return EXIT_SUCCESS;
}
