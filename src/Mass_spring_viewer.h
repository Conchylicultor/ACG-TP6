//=============================================================================
//
//   Exercise code for the lecture
//   "Advanced Computer Graphics"
//
//   Adapted from Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2013 LGG, epfl
//
//   DO NOT REDISTRIBUTE
//=============================================================================

#ifndef MASS_SPRING_VIEWER_H
#define MASS_SPRING_VIEWER_H


//== INCLUDES =================================================================

#include "Mass_spring_system.h"
#include "utils/Viewer_2D.h"
#include "utils/vec2.h"
#include <vector>
#include <fstream>

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>


//== CLASS DEFINITION =========================================================


class ImplicitSolver
{
public:
  ImplicitSolver ()
  { }

  void
  clear ()
  { triplets_.clear (); }

  void
  addElementToJacobian (const int row, const int col,
                        const float val)
  { triplets_.push_back (Eigen::Triplet<float> (row, col, val)); }

  void
  solve (float dt, float mass,
         std::vector<Particle> &particles);

private:
  Eigen::SparseLU<Eigen::SparseMatrix<float> > linear_solver_;
  std::vector<Eigen::Triplet<float> > triplets_;
};


/** \class Mass_spring_viewer Mass_spring_viewer.h <01-mass_springs/Mass_spring_viewer.h>
 Viewer class for the mass spring exercise.
 */
class Mass_spring_viewer : public Viewer_2D
{
public:
    /// constructor
    Mass_spring_viewer(const char* _title, int _width, int _height);
    /// destructor
    ~Mass_spring_viewer();

private: // GUI functions
    /// draw scene
    virtual void draw();

    /// handle keyboard events
    virtual void keyboard(int key, int x, int y);

    /// handle mouse events (used for interactive spring)
    virtual void mouse(int button, int state, int x, int y);

    /// handle mouse move events (used for interactive spring)
    virtual void motion(int x, int y);


private: // simulation functions
    /// compute all external and internal forces
    void compute_forces();

    /// perform one time step using either Euler, Midpoint, or Verlet
    void time_integration(float dt);

    /// perform impulse-based collision handling
    void impulse_based_collisions();

    void compute_jacobians(float dt);

private: // parameter settings
    /// parameter: mass of a particle
    float particle_mass_;

    /// parameter: amount of damping
    float damping_;

    /// parameter: strength of collision forces
    float collision_stiffness_;

    /// parameter: stiffness of springs
    float spring_stiffness_;

    /// parameter: internal damping of springs
    float spring_damping_;

    /// parameter: strength of area-preserving forces
    float area_stiffness_;

    /// paramters: which time-integration to use
    enum { Euler, Midpoint, Verlet, Implicit } integration_;

    /// paramters: external force type
    enum { None, Center, Gravitation,  } external_force_;

    /// parameter: how to handle collisiont
    enum { Force_based, Impulse_based } collisions_;

    /// parameter: use area-preserving forces?
    bool area_forces_;

    /// parameter: radius of particles (for rendering and collisions)
    float particle_radius_;

    /// parameter: visualize particle forces?
    bool show_forces_;

    /// parameter: visualize particle forces?
    bool equilibrium_forces_;

private: // simulation data
    /// the mass spring system to be simulated
    Mass_spring_system body_;

    /// class for storing the mouse spring
    struct Mouse_spring
    {
        /// position of the mouse cursor (one endpoint of spring)
        vec2 mouse_position;
        /// which particle is the other endpoint
        int  particle_index;
        /// is the spring active?
        bool active;
    };

    /// the interactive spring controlled by the mouse
    Mouse_spring mouse_spring_;

    ImplicitSolver solver_;

    // Collisions planes
    static const int NB_PANES = 4;
    static float planes[NB_PANES][3]; // Dirty C syntax
    vec2 planesNorms[NB_PANES]; // Norms of the plane
    float planesP[NB_PANES]; // Distance from the origin

    // For recording
    std::ofstream outputFile;
};


//=============================================================================
#endif // PHYSICS_ENGINE_VIEWER_HH defined
//=============================================================================

