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


//== INCLUDES =================================================================

#include "Mass_spring_viewer.h"
#include "utils/gl_text.h"
#include <sstream>


//== IMPLEMENTATION ==========================================================

float Mass_spring_viewer::planes[NB_PANES][3] = {
        {  0.0,  1.0, 1.0 },
        {  0.0, -1.0, 1.0 },
        {  1.0,  0.0, 1.0 },
        { -1.0,  0.0, 1.0 }/*,
        { -0.2, 1.0, 0.8 }*/ // Test with an oblique pane (Warning: change "NB_PANES" to 5);
    };

Mass_spring_viewer::
Mass_spring_viewer(const char* _title, int _width, int _height)
: Viewer_2D(_title, _width, _height)
{
    integration_         = Euler;
    collisions_          = Force_based;
    external_force_      = None;
    animate_             = false;
    area_forces_         = false;
    show_forces_         = false;
    equilibrium_forces_  = false;

    time_step_           = 0.001;
    particle_radius_     = 0.03;

    particle_mass_       = 0.1;
    damping_             = 0.1;
    collision_stiffness_ = 1000.0;
    spring_stiffness_    = 1000.0;
    spring_damping_      = 1.0;
    area_stiffness_      = 100000.0;

    mouse_spring_.active = false;


    // Collisions

    // Compute the panes normal & cie for the collisions
    for(unsigned int i=0; i<NB_PANES; ++i)
    {
        // According to http://mathworld.wolfram.com/HessianNormalForm.html

        float sqrtCoef = std::sqrt(planes[i][0]*planes[i][0] +
                                   planes[i][1]*planes[i][1]);

        planesNorms[i] = vec2(planes[i][0]/sqrtCoef,
                              planes[i][1]/sqrtCoef);

        planesP[i] = planes[i][2]/sqrtCoef;
    }
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::keyboard(int key, int x, int y)
{
    switch (key)
    {
        // setup problem 1
        case '1':
        {
            body_.clear();
            body_.add_particle( vec2(-0.5, -0.5), vec2(14.0, -2.0), particle_mass_, false );
            glutPostRedisplay();
            break;
        }

        // setup problem 2
        case '2':
        {
            body_.clear();
            for (int i=0; i<100; ++i)
            {
                body_.add_particle( vec2(0.9* cos(i/50.0*M_PI), 0.9*sin(i/50.0*M_PI)), vec2(-sin(i/50.0*M_PI), cos(i/50.0*M_PI)), particle_mass_, false );
            }

            glutPostRedisplay();
            break;
        }

        // setup problem 3
        case '3':
        {
            body_.clear();

            for (int i=0; i<10; ++i)
            {
                body_.add_particle( vec2(i*0.1, 0.8), vec2(0.0, 0.0), particle_mass_, i==0 );
            }

            for (unsigned int i=0; i<9; ++i)
            {
                body_.add_spring(i, i+1);
            }

            glutPostRedisplay();
            break;
        }

        // setup problem 4
        case '4':
        {
            body_.clear();

            body_.add_particle( vec2(-0.1, 0.7), vec2(0.0, 0.0), particle_mass_, false );
            body_.add_particle( vec2( 0.0, 0.6), vec2(0.0, 0.0), particle_mass_, false );
            body_.add_particle( vec2( 0.1, 0.7), vec2(0.0, 0.0), particle_mass_, false );

            body_.add_spring(0, 1);
            body_.add_spring(0, 2);
            body_.add_spring(1, 2);

            body_.add_triangle(0, 1, 2);

            glutPostRedisplay();
            break;
        }

        // setup problem 5
        case '5':
        {
            body_.clear();

            for (int i=0; i<8; ++i)
            {
                body_.add_particle( vec2(-0.5+0.2*cos(0.25*i*M_PI), -0.5+0.2*sin(0.25*i*M_PI)), vec2(5.0, 5.0), particle_mass_, false );
            }

            body_.add_particle( vec2(-0.5, -0.5), vec2(5.0, 5.0), particle_mass_, false );

            for (unsigned int i=0; i<8; ++i)
            {
                body_.add_spring( i, (i+1) % 8 );
                body_.add_spring(i, 8);
                body_.add_triangle(i, (i+1)%8, 8);
            }

            glutPostRedisplay();
            break;
        }

        // setup problem 5
        case '6':
        {
            body_.clear();

            body_.add_particle( vec2(0.0, 0.0), vec2(0.0, 0.0), particle_mass_, false );

            glutPostRedisplay();
            break;
        }

        // switch between time integration methods
        case 'i':
        {
            switch (integration_)
            {
                case Euler:    integration_ = Midpoint; break;
                case Midpoint: integration_ = Verlet;   break;
                case Verlet:   integration_ = Implicit;    break;
                case Implicit: integration_ = Euler; break;
            }
            glutPostRedisplay();
            break;
        }

        // switch between center and gravitation force
        case 'f':
        {
            switch (external_force_)
            {
                case None:        external_force_ = Center; break;
                case Center:      external_force_ = Gravitation; break;
                case Gravitation: external_force_ = None;      break;
            }
            glutPostRedisplay();
            break;
        }

        // switch between force-based and impulse-based collisions
        case 'c':
        {
            switch (collisions_)
            {
                case Force_based:   collisions_ = Impulse_based; break;
                case Impulse_based: collisions_ = Force_based;   break;
            }
            glutPostRedisplay();
            break;
        }


        // toggle area forces on/off
        case 'a':
        {
            area_forces_ = !area_forces_;
            glutPostRedisplay();
            break;
        }


        // visualization of particle forces on/off
        case 'v':
        {
            show_forces_ = !show_forces_;
            glutPostRedisplay();
            break;
        }


        // visualization of particle forces on/off
        case 'e':
        {
            equilibrium_forces_ = !equilibrium_forces_;
            glutPostRedisplay();
            break;
        }


        // let parent class do the work
        default:
        {
            Viewer_2D::keyboard(key, x, y);
            break;
        }
    }
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::draw()
{
    // parent's status text
    Viewer_2D::draw();


    // draw some status text
    glDisable(GL_LIGHTING);
    glColor3f(1,0,0);
    std::ostringstream oss;

    oss.str("");
    oss << "Integration: ";
    switch (integration_)
    {
        case Euler:    oss << "Euler";    break;
        case Midpoint: oss << "Midpoint"; break;
        case Verlet:   oss << "Verlet";   break;
        case Implicit: oss << "Implicit"; break;
    }
    glText(20, height_-40, oss.str());

    oss.str("");
    oss << "#Particles: " << body_.particles.size();
    glText(20, height_-60, oss.str());

    oss.str("");
    oss << "#Springs: " << body_.springs.size();
    glText(20, height_-80, oss.str());

    oss.str("");
    oss << "#Triangles: " << body_.triangles.size();
    glText(20, height_-100, oss.str());

    oss.str("");
    oss << "Area Forces: " << (area_forces_ ? "on" : "off");
    glText(20, height_-120, oss.str());

    oss.str("");
    oss << "Collisions: " << (collisions_ == Force_based ? "force" : "impulse");
    glText(20, height_-140, oss.str());

    oss.str("");
    oss << "External force: ";
    switch (external_force_)
    {
        case None:        oss << "None";        break;
        case Center:      oss << "Center";      break;
        case Gravitation: oss << "Gravitation"; break;
    }
    glText(20, height_-160, oss.str());

    oss.str("");
    oss << "Visualize forces: " << (show_forces_ ? "on" : "off");
    glText(20, height_-180, oss.str());

    oss.str("");
    oss << "Equilibrium forces: " << (equilibrium_forces_ ? "on" : "off");
    glText(20, height_-200, oss.str());


    // draw walls
    glDisable(GL_LIGHTING);
    glLineWidth(1.0);
    glColor3f(0.5,0.5,0.5);
    glBegin(GL_LINE_STRIP);
    glVertex2f( -1.0,  1.0 );
    glVertex2f( -1.0, -1.0 );
    glVertex2f(  1.0, -1.0 );
    glVertex2f(  1.0,  1.0 );
    glVertex2f( -1.0,  1.0 );
    glEnd();


    // draw mouse spring
    if (mouse_spring_.active)
    {
        glDisable(GL_LIGHTING);
        glLineWidth(5.0);
        glColor3f(1,0,0);
        glBegin(GL_LINES);
        glVertex2fv( body_.particles[mouse_spring_.particle_index].position.data() );
        glVertex2fv( mouse_spring_.mouse_position.data() );
        glEnd();
    }


    // draw particles, springs, triangles
    body_.draw(particle_radius_, show_forces_);
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::mouse(int _button, int _state, int _x, int _y)
{
    // need particles to do interaction
    if (!body_.particles.empty())
    {
        // mouse button release destroys current mouse spring
        if (_state == GLUT_UP)
        {
            mouse_spring_.active = false;
        }

        // mouse button press generates new mouse spring
        else if (_state == GLUT_DOWN)
        {
            // get point under mouse cursor
            vec2 p = pick(_x, _y);

            // find closest particle
            int   pidx = -1;
            float dmin = FLT_MAX;
            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                float d = norm(p - body_.particles[i].position);
                if (d < dmin)
                {
                    dmin = d;
                    pidx = i;
                }
            }

            // construct mouse spring
            mouse_spring_.mouse_position = p;
            mouse_spring_.particle_index = pidx;
            mouse_spring_.active = true;
        }
    }

    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::motion(int _x, int _y)
{
    if (mouse_spring_.active)
    {
        // update mouse positions
        mouse_spring_.mouse_position = pick(_x, _y);
    }

    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::time_integration(float dt)
{
    switch (integration_)
    {
        case Euler:
        {
            /** \todo (Part 1) Implement Euler integration scheme
             \li The Particle class has variables position_t and velocity_t to store current values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */

            // Compute forces
            compute_forces();

            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                if(body_.particles[i].locked)
                {
                    continue; // No change
                }
                body_.particles[i].acceleration = body_.particles[i].force / body_.particles[i].mass;
                body_.particles[i].velocity += dt*body_.particles[i].acceleration;
                body_.particles[i].position += dt*body_.particles[i].velocity; // Use new velocity
            }

            break;
        }

        case Midpoint:
        {
            /** \todo (Part 2) Implement the Midpoint time integration scheme
             \li The Particle class has variables position_t and velocity_t to store current values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */

            compute_forces(); // For x1, v1

            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                if(body_.particles[i].locked)
                {
                    continue; // No change
                }

                // Save current pos and speed
                body_.particles[i].position_t = body_.particles[i].position; // Save x(t)
                body_.particles[i].velocity_t = body_.particles[i].velocity; // Save v(t)

                // Compute ...
                body_.particles[i].acceleration = body_.particles[i].force / body_.particles[i].mass; // a1
                body_.particles[i].position += dt/2*body_.particles[i].velocity; // x2
                body_.particles[i].velocity += dt/2*body_.particles[i].acceleration; // v2

            }

            compute_forces(); // For x2, v2

            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                if(body_.particles[i].locked)
                {
                    continue; // No change
                }

                body_.particles[i].acceleration = body_.particles[i].force / body_.particles[i].mass; // a2
                body_.particles[i].position = body_.particles[i].position_t + dt/2*body_.particles[i].velocity; // x(t+h)
                body_.particles[i].velocity = body_.particles[i].velocity_t + dt/2*body_.particles[i].acceleration; // v(t+h)
            }

            break;
        }


        case Verlet:
        {
            /** \todo (Part 2) Implement the Verlet time integration scheme
             \li The Particle class has a variable acceleration to remember the previous values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */

            compute_forces();

            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                if(body_.particles[i].locked)
                {
                    continue; // No change
                }

                body_.particles[i].acceleration = body_.particles[i].force / body_.particles[i].mass; // a(t)
                body_.particles[i].position += dt*body_.particles[i].velocity + ((dt*dt)/2)*body_.particles[i].acceleration; // x(t+h)
            }

            compute_forces();

            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                if(body_.particles[i].locked)
                {
                    continue; // No change
                }

                vec2 aTPlusH = body_.particles[i].force / body_.particles[i].mass; // a(t+h)
                body_.particles[i].velocity += dt*(body_.particles[i].acceleration + aTPlusH) / 2; // v(t+h)
            }

            break;
        }

        case Implicit:
        {
            /// The usual force computation method is called, and then the jacobian matrix dF/dx is calculated
            compute_forces ();
            compute_jacobians (dt);

            /// Finally the linear system is composed and solved
            solver_.solve (dt, particle_mass_, body_.particles);

            break;
        }
    }


    // impulse-based collision handling
    if (collisions_ == Impulse_based)
    {
        impulse_based_collisions();
    }


    glutPostRedisplay();
}


//-----------------------------------------------------------------------------

// Area of the triangle
float computeAreaTriangle(const Triangle &triangle)
{
    const vec2 &pt1 = triangle.particle0->position;
    const vec2 &pt2 = triangle.particle1->position;
    const vec2 &pt3 = triangle.particle2->position;
    // Using the Shoelace formula:
    return 1.0f/2.0f * std::abs((pt1[0] - pt3[0])*(pt2[1] - pt1[1]) -
                                (pt1[0] - pt2[0])*(pt3[1] - pt1[1]));
}

// Partial derivate with respect to pt1
vec2 computeDerivateAreaTriangle(const vec2 &pt1, const vec2 &pt2, const vec2 &pt3)
{
    // We differenciate the triangle area:
    /*
    float u = pt1[0]*pt2[1] + pt2[0]*pt3[1] + pt3[0]*pt1[1]
            - pt1[0]*pt3[1] - pt3[0]*pt2[1] - pt2[0]*pt1[1];
    float sign = u / std::abs(u);
    */
    float sign = 1.0;
    return 1.0f/2.0f * sign * vec2((pt2[1] - pt3[1]),
                                   (pt3[0] - pt2[0]));
}

void
Mass_spring_viewer::compute_forces()
{
    // clear forces
    for (unsigned int i=0; i<body_.particles.size(); ++i)
        body_.particles[i].force = vec2(0,0);


    /** \todo (Part 1) Implement center force
     */
    if (external_force_ == Center)
    {
        for (unsigned int i=0; i<body_.particles.size(); ++i)
        {
            body_.particles[i].force += 20.0f*(vec2(0,0) - body_.particles[i].position);
        }
    }


    /** \todo (Part 1) Implement damping force
     \li The damping coefficient is given as member damping_
     */
    /// Do damping only for explicit methods
    if (integration_ != Implicit)
      for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
        p_it->force -= damping_ * p_it->velocity;


    /** \todo (Part 1) Implement gravitation force
     \li Particle mass available as particle_mass_
     */
    if (external_force_ == Gravitation)
    {
        for (unsigned int i=0; i<body_.particles.size(); ++i)
        {
            body_.particles[i].force += vec2(0.f,-9.81)*body_.particles[i].mass;
            // body_.particles[i].force += vec2(0.f,-9.81)*particle_mass_;
        }
    }


    /** \todo (Part 1) Implement force based boundary collisions
     \li Collision coefficient given as collision_stiffness_
     */
    // collision forces
    if (collisions_ == Force_based)
    {
        for (unsigned int i=0; i<body_.particles.size(); ++i)
        {
            for (unsigned int j=0 ; j<NB_PANES ; ++j)
            {
                // From:
                // - Pane equation: ax + by + c = 0
                // - Point (x0,y0)
                // We compute the distance according to the formula:
                // - Dist = (ax0 + by0 + c) / sqrt(a^2 + b^2)

                float distPointPlanes = dot(planesNorms[j], body_.particles[i].position) + planesP[j];
                if(distPointPlanes < particle_radius_) // Collision
                {
                    body_.particles[i].force += collision_stiffness_ * (particle_radius_-distPointPlanes) * planesNorms[j];
                }
            }
        }
    }


    /** \todo (Part 1) Compute force of the interactive mass spring
     \li Required coefficients are given as spring_stiffness_ and spring_damping_
     */
    if (mouse_spring_.active)
    {
        Particle& p0 = body_.particles[ mouse_spring_.particle_index ];

        vec2 pos0 = p0.position;
        vec2 pos1 = mouse_spring_.mouse_position;

        vec2 vel0 = p0.velocity;
        vec2 vel1 = vec2(0,0);

        float d = norm(pos0 - pos1);

        p0.force += -1 * (spring_stiffness_*d + spring_damping_*dot(vel0-vel1, pos0-pos1)/d) * ((pos0-pos1)/d);

    }


    /** \todo (Part 1) Compute spring forces
     \li Required information about springs in the scene are found in body_.springs
     \li Required coefficients are given as spring_stiffness_ and spring_damping_
     */
    for (unsigned int i=0; i<body_.springs.size(); ++i){
        Spring &nextSpring = body_.springs[i];

        vec2 pos0 = nextSpring.particle0->position;
        vec2 pos1 = nextSpring.particle1->position;

        vec2 vel0 = nextSpring.particle0->velocity;
        vec2 vel1 = nextSpring.particle1->velocity;

        float d = norm(pos0 - pos1);

        vec2 f0 = -1 * (spring_stiffness_*(d-nextSpring.rest_length) + spring_damping_*dot(vel0-vel1, pos0-pos1)/d) * ((pos0-pos1)/d);
        nextSpring.particle0->force += f0;
        nextSpring.particle1->force -= f0;
    }


    /** \todo (Part 2) Compute more forces in part 2 of the exercise: triangle-area forces, binding forces, etc.
     */
    if (area_forces_)
    {
        for (unsigned int k = 0; k < body_.triangles.size(); k++) {
            float ka = area_stiffness_;
            Triangle &triangle = body_.triangles[k];

            // float area = computeAreaTriangle(triangle); // Right but useless
            // std::cout << triangle.area() << " " << area << std::endl; // Because

            float A = triangle.rest_area;
            float EaFixed = ka * (triangle.area() - A);

            // Forces
            const vec2 &pt1 = triangle.particle0->position;
            const vec2 &pt2 = triangle.particle1->position;
            const vec2 &pt3 = triangle.particle2->position;

            triangle.particle0->force += - EaFixed * computeDerivateAreaTriangle(pt1, pt2, pt3); // The 1/2 factor it inside the computeDerivateAreaTriangle
            triangle.particle1->force += - EaFixed * computeDerivateAreaTriangle(pt2, pt3, pt1);
            triangle.particle2->force += - EaFixed * computeDerivateAreaTriangle(pt3, pt1, pt2);
        }
    }

    // Equilibrium forces
    if (equilibrium_forces_ && integration_ != Implicit)
    {
        float coeffRepulsion = 0.02f;
        for (unsigned int i = 0; i < body_.particles.size(); i++)
        {
            for (unsigned int j = 0; j < body_.particles.size(); j++)
            {
                if(i != j)
                {
                    const vec2 &pt1 = body_.particles.at(i).position;
                    const vec2 &pt2 = body_.particles.at(j).position;

                    float dist = norm(pt2 - pt1);
                    vec2 direction = (pt2-pt1)/dist;

                    body_.particles.at(i).force -= coeffRepulsion / (dist*dist) * direction;
                }
            }
        }

    }
}

//-----------------------------------------------------------------------------


void Mass_spring_viewer::impulse_based_collisions()
{
    /** \todo (Part 2) Handle collisions based on impulses
     */

    const float restitutionCoef = 1.0;
    for (unsigned int i=0; i<body_.particles.size(); ++i)
    {
        // Detect collision

        for (unsigned int j=0 ; j<NB_PANES ; ++j)
        {
            // From:
            // - The pane equation: ax + by + c = 0
            // - Our particule position: Point (x0,y0)
            // We compute the distance according to the formula:
            // - DistPointPane = (a*x0 + b*y0 + c) / sqrt(a^2 + b^2)

            float distPointPlanes = dot(planesNorms[j], body_.particles[i].position) + planesP[j];
            if(distPointPlanes < particle_radius_ && dot(body_.particles[i].velocity,planesNorms[j]) < 0) // Collision && check direction
            {
                // j = -(1+e)p.n
                float jCoeff = -(1+restitutionCoef)*dot(body_.particles[i].velocity, planesNorms[j]);
                // v(t+dt) = v(t) + j*n
                body_.particles[i].velocity += jCoeff*planesNorms[j];
            }
        }
    }
}

//=============================================================================

void computeJacobian_TriangleArea (const vec2 &Fpt0, const vec2 &pt1, const vec2 &pt2, int row, int col,float C, ImplicitSolver &solver_){
    // compute the Jacobian of triangle area for one force Fpt0 and one particule pt0

    float DFxDptx = Fpt0[0]*(pt1[1]-pt2[1])/C;
    float DFxDpty = Fpt0[0]*(pt2[0]-pt1[0])/C;
    float DFyDptx = Fpt0[1]*(pt1[1]-pt2[1])/C;
    float DFyDpty = Fpt0[1]*(pt2[0]-pt1[0])/C;

    solver_.addElementToJacobian(row,col,DFxDptx);
    solver_.addElementToJacobian(row,col+1,DFxDpty);
    solver_.addElementToJacobian(row+1,col,DFyDptx);
    solver_.addElementToJacobian(row+1,col+1,DFyDpty);
}

//=============================================================================

void Mass_spring_viewer::compute_jacobians (float dt)
{
    /// Clear the solver matrices
    solver_.clear ();

    /** \todo (Part 2) Implement the corresponding jacobians for each of the force types.
   * Use the code from compute_forces() as the starting ground.
   */

    // Mouse spring
    if (mouse_spring_.active)
    {
        // http://blog.mmacklin.com/2012/05/04/implicitsprings/

        Particle& part_i = body_.particles[mouse_spring_.particle_index];

        vec2 pij = part_i.position - mouse_spring_.mouse_position;
        if(norm(pij) != 0) // Avoid singularity
        {
            vec2 pijn = pij / norm(pij);

            std::cout << pij << std::endl;
            std::cout << pijn << std::endl;

            Eigen::Vector2f xij(pijn[0], pijn[1]);
            Eigen::Matrix2f xijxijt = xij*xij.transpose();
            //Eigen::Matrix2f dFdxi = -spring_stiffness_ * ( (1-0/norm(pij)) * (Eigen::Matrix2f::Identity()-xijxijt) + xijxijt ); // Spring term
            Eigen::Matrix2f dFdxi = -spring_stiffness_ * Eigen::Matrix2f::Identity(); // Spring term
            dFdxi                += -spring_damping_   * 1/dt * xijxijt;// Dampling term

            std::cout << xij << std::endl;
            std::cout << xijxijt << std::endl;
            std::cout << dFdxi << std::endl;

            unsigned int y = 2 * mouse_spring_.particle_index;
            solver_.addElementToJacobian(y,   y,   dFdxi(0,0));
            solver_.addElementToJacobian(y,   y+1, dFdxi(0,1));
            solver_.addElementToJacobian(y+1, y  , dFdxi(1,0));
            solver_.addElementToJacobian(y+1, y+1, dFdxi(1,1));
        }
    }

    // Damped springs
    for (unsigned int i=0; i<body_.springs.size(); ++i){

        Spring &nextSpring = body_.springs[i];
        Particle *part_i = nextSpring.particle0;
        Particle *part_j = nextSpring.particle1;

        vec2 pij = part_i->position - part_j->position;
        vec2 pijn = pij / norm(pij);

        Eigen::Vector2f xij(pijn[0], pijn[1]);
        Eigen::Matrix2f xijxijt = xij*xij.transpose();
        Eigen::Matrix2f dFdxi = -spring_stiffness_ * ( (1-nextSpring.rest_length/norm(pij)) * (Eigen::Matrix2f::Identity()-xijxijt) + xijxijt ); // Spring term
        dFdxi                += -spring_damping_   * 1/dt * xijxijt; // Dampling term

        unsigned int yi = 2 * part_i->id;
        unsigned int yj = 2 * part_j->id;

        solver_.addElementToJacobian(yi,   yi,   dFdxi(0,0));
        solver_.addElementToJacobian(yi,   yi+1, dFdxi(0,1));
        solver_.addElementToJacobian(yi+1, yi  , dFdxi(1,0));
        solver_.addElementToJacobian(yi+1, yi+1, dFdxi(1,1));

        solver_.addElementToJacobian(yj,   yj,   dFdxi(0,0));
        solver_.addElementToJacobian(yj,   yj+1, dFdxi(0,1));
        solver_.addElementToJacobian(yj+1, yj  , dFdxi(1,0));
        solver_.addElementToJacobian(yj+1, yj+1, dFdxi(1,1));

        dFdxi = -dFdxi;

        solver_.addElementToJacobian(yi,   yj,   dFdxi(0,0));
        solver_.addElementToJacobian(yi,   yj+1, dFdxi(0,1));
        solver_.addElementToJacobian(yi+1, yj  , dFdxi(1,0));
        solver_.addElementToJacobian(yi+1, yj+1, dFdxi(1,1));

        solver_.addElementToJacobian(yj,   yi,   dFdxi(0,0));
        solver_.addElementToJacobian(yj,   yi+1, dFdxi(0,1));
        solver_.addElementToJacobian(yj+1, yi  , dFdxi(1,0));
        solver_.addElementToJacobian(yj+1, yi+1, dFdxi(1,1));
    }

    // Force based collisions
    if (collisions_ == Force_based)
    {
        for (unsigned int i=0; i<body_.particles.size(); ++i)
        {
            for (unsigned int j=0 ; j<NB_PANES ; ++j)
            {
                // From:
                // - Pane equation: ax + by + c = 0
                // - Point (x0,y0)
                // We compute the distance according to the formula:
                // - Dist = (ax0 + by0 + c) / sqrt(a^2 + b^2)

                float distPointPlanes = dot(planesNorms[j], body_.particles[i].position) + planesP[j];
                if(distPointPlanes < particle_radius_) // Collision
                {
                    vec2 dFdx = - collision_stiffness_ * planesNorms[j][0] * planesNorms[j];
                    vec2 dFdy = - collision_stiffness_ * planesNorms[j][1] * planesNorms[j];

                    solver_.addElementToJacobian(2*i,  2*i,   dFdx[0]);
                    solver_.addElementToJacobian(2*i,  2*i+1, dFdy[0]);
                    solver_.addElementToJacobian(2*i+1,2*i,   dFdx[1]);
                    solver_.addElementToJacobian(2*i+1,2*i+1, dFdy[1]);
                }
            }
        }
    }

    // Gravitation
    // Constant force so jacobian = 0

    // Center force
    if (external_force_ == Center)
    {
        const float centerForceStrength = 20.0f;
        for (unsigned int i=0; i<body_.particles.size(); ++i)
        {
            solver_.addElementToJacobian(2*i,2*i, -1.0f * centerForceStrength);
            solver_.addElementToJacobian(2*i+1,2*i+1, -1.0f * centerForceStrength);
        }
    }

    // Triangle area forces
    if (area_forces_)
    {
        //For each triangle
        for (unsigned int k = 0; k < body_.triangles.size(); k++) {
            float ka = area_stiffness_;
            Triangle &triangle = body_.triangles[k];

            float C = triangle.area()-triangle.rest_area;
            // Sommets
            Particle* pt[3] = {triangle.particle0, triangle.particle1, triangle.particle2};


            vec2 Fpt[3];
            Fpt[0] = - ka*C * computeDerivateAreaTriangle(pt[0]->position, pt[1]->position, pt[2]->position);
            Fpt[1] = - ka*C * computeDerivateAreaTriangle(pt[1]->position, pt[2]->position, pt[0]->position);
            Fpt[2] = - ka*C * computeDerivateAreaTriangle(pt[2]->position, pt[0]->position, pt[1]->position);

            // Double boucle : pour chaque force Fj et pour chaque particule pi
            for (int pi = 0 ; pi<3 ; pi++) {
                for (int Fj = 0 ; Fj<3; Fj++) {
                    int row = 2*(pt[Fj]->id);
                    int col = 2*(pt[pi%3]->id);
                    computeJacobian_TriangleArea(Fpt[Fj],pt[(pi+1)%3]->position,pt[(pi+2)%3]->position,row,col,C,solver_);
                }
            }

        }
    }
}

//=============================================================================
void ImplicitSolver::solve (float dt, float mass,
                            std::vector<Particle> &particles)
{
  int num_particles = particles.size ();

  /// Build the Jacobian matrix from the sparse set of elements
  Eigen::SparseMatrix<float> J (2 * num_particles, 2 * num_particles);
  J.setFromTriplets (triplets_.begin (), triplets_.end ());

  /// Build up the position, velocity and force vectors
  Eigen::VectorXf pos_vec (Eigen::VectorXf::Zero (2 * num_particles)),
                  velocity_vec (Eigen::VectorXf::Zero (2 * num_particles)),
                  force_vec (Eigen::VectorXf::Zero (2 * num_particles));
  for (size_t p_i = 0; p_i < num_particles; ++p_i)
  {
    pos_vec (2 * p_i + 0) = particles[p_i].position[0];
    pos_vec (2 * p_i + 1) = particles[p_i].position[1];
    velocity_vec (2 * p_i + 0) = particles[p_i].velocity[0];
    velocity_vec (2 * p_i + 1) = particles[p_i].velocity[1];
    force_vec (2 * p_i + 0) = particles[p_i].force[0];
    force_vec (2 * p_i + 1) = particles[p_i].force[1];
  }

  /// Kick out the fixed particles by creating a sparse selection matrix
  std::vector<Eigen::Triplet<float> > triplets_selection;
  int valid_particle_index = 0;
  for (size_t p_i = 0; p_i < num_particles; ++p_i)
    if (!particles[p_i].locked)
    {
      triplets_selection.push_back (Eigen::Triplet<float> (2 * p_i + 0, 2 * valid_particle_index + 0, 1.f));
      triplets_selection.push_back (Eigen::Triplet<float> (2 * p_i + 1, 2 * valid_particle_index + 1, 1.f));
      valid_particle_index ++;
    }
  Eigen::SparseMatrix<float> mat_selection (2 * num_particles, 2 * valid_particle_index);
  mat_selection.setFromTriplets (triplets_selection.begin (), triplets_selection.end ());

  /// Sparse identity matrix
  Eigen::SparseMatrix<float> Id (2 * valid_particle_index, 2 * valid_particle_index);
  Id.setIdentity ();

  /// Apply the selection matrix on each vector and the Jacobian matrix
  pos_vec = mat_selection.transpose () * pos_vec;
  velocity_vec = mat_selection.transpose () * velocity_vec;
  force_vec = mat_selection.transpose () * force_vec;
  J = mat_selection.transpose () * J * mat_selection;

  /// Build the right and left hand sides of the linear system
  Eigen::SparseMatrix<float> A = Id - dt * dt / mass * J;
  Eigen::VectorXf b;
  b = dt * velocity_vec + dt * dt / mass * force_vec + (Id - dt * dt / mass * J) * pos_vec;

  /// Solve the system and use the selection matrix again to arrange the new positions in a vector
  linear_solver_.analyzePattern (A);
  linear_solver_.compute (A);
  Eigen::VectorXf new_pos = mat_selection * linear_solver_.solve (b);

  /// Extract the positions from the solution vector and set the new positions and velocities inside the particle structures
  for (size_t p_i = 0; p_i < num_particles; ++p_i)
  {
    if (!particles[p_i].locked)
    {
      vec2 pos_update (new_pos (2 * p_i + 0), new_pos (2 * p_i + 1));
      particles[p_i].velocity = (pos_update - particles[p_i].position) / dt;
      particles[p_i].position = pos_update;
    }
  }
}
