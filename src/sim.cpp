#include <vector>
#include <unordered_map>
#include <string>
#include <math.h>
#include "particle.cpp"

using namespace std;

#define PI 3.14159265 // Should be used from mathlib

class Force {
    public:
        float magnitude;
        Vector3d dir;
        Force(float mag, Vector3d dir1): magnitude(mag), dir(dir1) {}
};
class Simulator {
    public:
        class VectorHasher {
            public:
                Simulator* sim;
                VectorHasher(Simulator* s): sim(s) {}
                size_t operator()(Vector3d v) const {
                    int x = (int) (v.x / sim->d);
                    int y = (int) (v.y / sim->d);
                    int z = (int) (v.z / sim->d);

                    return ((hash<int>()(x)
                    ^ (hash<int>()(y) << 1)) >> 1)
                    ^ (hash<int>()(z) << 1);
                }
        };
        class VectorEqual {
            public:
                Simulator* sim;
                VectorEqual(Simulator* s): sim(s) {}
                size_t operator()(Vector3d v1, Vector3d v2) const {
                    int x1 = (int) (v1.x / sim->d);
                    int y1 = (int) (v1.y / sim->d);
                    int z1 = (int) (v1.z / sim->d);

                    int x2 = (int) (v2.x / sim->d);
                    int y2 = (int) (v2.y / sim->d);
                    int z2 = (int) (v2.z / sim->d);

                    return (x1 == x2) && (y1 == y2) && (z1 == z2);
                }
        };
        vector<Particle*> particles;
        float g = 0.98f * 2;
        float r0 = 1.0f;
        float k = 8.314;
        float gamma = 1;
        float mu = 0.5;
        float xl = 0;
        float yl = 0;
        float zl = 0;
        float xt;
        float yt;
        float zt;
        unordered_map<Vector3d,vector<Particle*>*, VectorHasher, VectorEqual>* cells;
        Simulator(int x, int y, int z, int xbounds, int ybounds, int zbounds, float din=1) {
            for (int i = 0; i < x; i++) {
                for (int j = 0; j < y; j++) {
                    for (int k = z/2; k < z; k++) {
                        Vector3d pos(xbounds * 1.0/(x + 1) * (i + 1), ybounds * 1.0/(y + 1) * (j + 1), 
                            zbounds * 1.0/(z + 1) * (k + 1));
                        Particle* p = new Particle(pos);
                        p->index = particles.size();
                        particles.push_back(p);
                    }
                }
            }
            d = din;
            poly6Const = 315/(64 * PI * pow(d, 9));
            spikyGradConst = -45/(PI * pow(d, 6));
            d2 = d * d;

            xt = xbounds;
            yt = ybounds;
            zt = zbounds;

            buildGhosts();

            cells = new unordered_map<Vector3d,vector<Particle*>*, VectorHasher, VectorEqual>(100, VectorHasher(this), VectorEqual(this));
        }
        float getDensity(Vector3d pos, bool ghost=false) {
            vector<Particle*> neighborhood = getNeighborhood(pos);
            float density = 0;
            for (int i = 0; i < neighborhood.size(); i++) {
                if (neighborhood[i]->ghost == ghost) {
                    Vector3d dist = Vector3d::sub(neighborhood[i]->pos, pos);
                    density += neighborhood[i]->mass * poly6(dist.magnitude());
                }
            }
            return density;
        }
        void stepSPH(float timestep=0.01f) {
            vector<Force*> gravityF(particles.size());
            vector<Force*> viscosityF(particles.size());
            vector<Force*> pressureF(particles.size());

            vector<vector<Particle*>> neighborhoods(particles.size());
            //Calculate neighborhoods
            updateNeighborhoods();
            for (int i = 0; i < particles.size(); i++) {
                neighborhoods[i] = getNeighborhood(particles[i]->pos);
            }
            //Calculate densities and pressures
            vector<float> density(particles.size());
            vector<float> pressure(particles.size());
            for (int i = 0; i < particles.size(); i++) {
                for (int j = 0; j < neighborhoods[i].size(); j++) {
                    Vector3d dist = Vector3d::sub(neighborhoods[i][j]->pos, particles[i]->pos);
                    density[i] += neighborhoods[i][j]->mass * poly6(dist.magnitude());
                }   

                //calculate pressure
                pressure[i] = k * r0/gamma * (pow(density[i]/r0, gamma) - 1);
            }
            //Calculate pressure force
            for (int i = 0; i < particles.size(); i++) {
                Vector3d pressureT;
                for (int j = 0; j < neighborhoods[i].size(); j++) {
                    
                    Vector3d dist = Vector3d::sub(neighborhoods[i][j]->pos, particles[i]->pos);
                    float pIndex = neighborhoods[i][j]->index;
                    float constant = neighborhoods[i][j]->mass * (pressure[i] + pressure[pIndex])/(2 * density[pIndex]);
                    if (density[pIndex] == 0) {
                        constant = 0;
                    }
                    /*printf("Density: %f\n", density[pIndex]);
                    printf("pressure: %f\n", pressure[pIndex]);
                    printf("pressureme: %f\n", pressure[i]);
                    printf("Constant: %f\n", constant);*/
                    pressureT = Vector3d::add(Vector3d::mult_num(spikyGrad(dist), constant), pressureT);
                }
                pressureT = Vector3d::mult_num(pressureT, 1);
                pressureF[i] = new Force(1, pressureT);
            }
            //Calculate viscosity force
            for (int i = 0; i < particles.size(); i++) {
                Vector3d viscosityT;
                for (int j = 0; j < neighborhoods[i].size(); j++) {
                    Vector3d dist = Vector3d::sub(neighborhoods[i][j]->pos, particles[i]->pos);
                    float pIndex = neighborhoods[i][j]->index;
                    float constant = neighborhoods[i][j]->mass * laplacian(dist.magnitude()) / density[pIndex];
                    if (density[pIndex] == 0) {
                        constant = 0;
                    }
                    Vector3d velDiff = Vector3d::sub(neighborhoods[i][j]->vel, particles[i]->vel);
                    viscosityT = Vector3d::add(Vector3d::mult_num(velDiff, constant), viscosityT);
                }

                viscosityT = Vector3d::mult_num(viscosityT, mu);
                viscosityF[i] = new Force(1, viscosityT);
            }
            //Calculate gravitational effects
            for (int i = 0; i < particles.size(); i++) {
                float mag = density[i] * g;
                gravityF[i] = new Force(mag, Vector3d(0, 0, -1));
            }
            for (int i = 0; i < particles.size(); i++) {
                if (particles[i]->ghost == false) {
                    Vector3d totalForce = Vector3d::mult_num(gravityF[i]->dir, gravityF[i]->magnitude);
                    Vector3d pressureVec = Vector3d::mult_num(pressureF[i]->dir, pressureF[i]->magnitude);
                    Vector3d viscosityVec = Vector3d::mult_num(viscosityF[i]->dir, viscosityF[i]->magnitude);
                    totalForce = Vector3d::add(totalForce, pressureVec);
                    totalForce = Vector3d::add(totalForce, viscosityVec);

                    //Calculate acceleration
                    if (density[i] != 0) {
                        Vector3d acceleration = Vector3d::mult_num(totalForce, 1.0/density[i]);
                        particles[i]->updateVel(acceleration, timestep);
                        particles[i]->updatePos(particles[i]->vel, timestep);
                        fixPosition(particles[i]);
                    }
                }
            }
        }
    private:
        float d;
        float poly6Const;
        float spikyGradConst;
        float d2;
        float poly6(float r) {
            if (r > d || r < 0) {
                return 0;
            }

            return poly6Const * pow((d2 - r*r), 3.0);
        }
        Vector3d spikyGrad(Vector3d r) {
            float mag = r.magnitude();
            if (mag == 0) {
                return Vector3d();
            }
            if (mag > d || mag < 0) {
                return Vector3d();
            }
            r.normalize();
            return Vector3d::mult_num(r, spikyGradConst * pow((d - mag), 2));
        }
        float laplacian(float r) {
            if (r > d || r < 0) {
                return 0;
            }
            return -spikyGradConst * (d - r);
        }
        void buildGhosts() {
            //build cube
            float increment = 0.2;
            //x == 0
            for (float i = 0; i < yt; i+=increment) {
                for (float j = 0; j < zt; j+=increment) {
                    Vector3d pos(0, i, j);
                    Particle* p = new Particle(pos);
                    p->ghost = true;
                    p->index = particles.size();
                    particles.push_back(p);
                }
            }
            //x == xt
            for (float i = 0; i < yt; i+=increment) {
                for (float j = 0; j < zt; j+=increment) {
                    Vector3d pos(xt, i, j);
                    Particle* p = new Particle(pos);
                    p->ghost = true;
                    p->index = particles.size();
                    particles.push_back(p);
                }
            }

            //y = 0
            /*for (float i = 0; i < xt; i+=increment) {
                for (float j = 0; j < zt; j+=increment) {
                    Vector3d pos(i, 0, j);
                    Particle* p = new Particle(pos);
                    p->ghost = true;
                    p->index = particles.size();
                    particles.push_back(p);
                }
            }
            //y == yt
            for (float i = 0; i < xt; i+=increment) {
                for (float j = 0; j < zt; j+=increment) {
                    Vector3d pos(i, yt, j);
                    Particle* p = new Particle(pos);
                    p->ghost = true;
                    p->index = particles.size();
                    particles.push_back(p);
                }
            }*/
            
            //z = 0
            for (float i = 0; i < xt; i+=increment) {
                for (float j = 0; j < yt; j+=increment) {
                    Vector3d pos(i, j, 0);
                    Particle* p = new Particle(pos);
                    p->ghost = true;
                    p->index = particles.size();
                    particles.push_back(p);
                }
            }
            //z = zt
            for (float i = 0; i < xt; i+=increment) {
                for (float j = 0; j < yt; j+=increment) {
                    Vector3d pos(i, j, zt);
                    Particle* p = new Particle(pos);
                    p->ghost = true;
                    p->index = particles.size();
                    particles.push_back(p);
                }
            }
        }
        void updateNeighborhoods() {
            unordered_map<Vector3d,vector<Particle*>*, VectorHasher, VectorEqual>::iterator it;
            for (it = cells->begin(); it != cells->end(); it++) {
                delete it->second;
            }
            cells -> clear();
            for (int i = 0; i < particles.size(); i++) {
                Vector3d location = particles[i]->pos;
                if (cells->count(location) == 0) {
                    cells->emplace(location, new vector<Particle*>());
                }
                cells->at(location)->push_back(particles[i]);
            }
        }
        void updateNeighbors(Vector3d location, vector<Particle*>* neighbors) {
            if (cells->count(location) != 0) {
                vector<Particle*>* neighborhood = cells->at(location);
                neighbors->insert(neighbors->end(), neighborhood->begin(), neighborhood->end());
            }
        }
        vector<Particle*> getNeighborhood(Vector3d location) {
            Vector3d xLM(location.x - d, location.y, location.z);
            Vector3d xRM(location.x + d, location.y, location.z);
            Vector3d yLM(location.x, location.y - d, location.z);
            Vector3d yRM(location.x, location.y + d, location.z);
            Vector3d zLM(location.x, location.y, location.z - d);
            Vector3d zRM(location.x, location.y, location.z + d);

            /*Vector3d xyL(location.x - d, location.y - d, location.z);
            Vector3d xyR(location.x + d, location.y - d, location.z);
            Vector3d yxL(location.x - d, location.y + d, location.z);
            Vector3d yxR(location.x + d, location.y + d, location.z);

            Vector3d xzL(location.x - d, location.y, location.z - d);
            Vector3d xzR(location.x + d, location.y, location.z - d);
            Vector3d zxL(location.x - d, location.y, location.z + d);
            Vector3d zxR(location.x + d, location.y, location.z + d);

            Vector3d yzL(location.x, location.y - d, location.z - d);
            Vector3d yzR(location.x, location.y + d, location.z - d);
            Vector3d zyL(location.x, location.y - d, location.z + d);
            Vector3d zyR(location.x, location.y + d, location.z + d);*/
            
            vector<Particle*> neighbors;
            updateNeighbors(location, &neighbors);
            updateNeighbors(xLM, &neighbors);
            updateNeighbors(xRM, &neighbors);
            updateNeighbors(yLM, &neighbors);
            updateNeighbors(yRM, &neighbors);
            updateNeighbors(zLM, &neighbors);
            updateNeighbors(zRM, &neighbors);
            return neighbors;//FIX THIS
        }
        void fixPosition(Particle* particle) {
            if (particle->pos.x > xt-0.1) {
                particle->pos.x = xt-0.1;
            }
            if (particle->pos.x < xl + 0.1) {
                particle->pos.x = xl + 0.1;
            }
            if (particle->pos.y > yt - 0.1) {
                particle->pos.y = yt - 0.1;
            }
            if (particle->pos.y < yl + 0.1) {
                particle->pos.y = yl + 0.1;
            }
            if (particle->pos.z > zt - 0.1) {
                particle->pos.z = zt - 0.1;
            }
            if (particle->pos.z < zl + 0.1) {
                particle->pos.z = zl + 0.1;
            }
        }
};
