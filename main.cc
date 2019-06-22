
#include <GL/freeglut.h>
#include <sys/sysinfo.h>
#include "Common.h"
#include "pcg_basic.h"
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <vector>

#include <memory>
#include <cstdlib>
#include <unistd.h>
#include <map>
#include <vector>
#include <climits>
#include <list>
#include <atomic>
#include <set>

#include <pthread.h>

using namespace std;

const int MaxLeaf = 8;

void DrawBox(Extrema2i const&area) {
  Vec2i pt;
  glBegin(GL_LINE_STRIP);
  pt = area.mMin + Vec2i(0,0) * area.GetSize();
  glVertex2iv(&pt[0]);
  pt = area.mMin + Vec2i(0,1) * area.GetSize();
  glVertex2iv(&pt[0]);
  pt = area.mMin + Vec2i(1,1) * area.GetSize();
  glVertex2iv(&pt[0]);
  pt = area.mMin + Vec2i(1,0) * area.GetSize();
  glVertex2iv(&pt[0]);
  pt = area.mMin + Vec2i(0,0) * area.GetSize();
  glVertex2iv(&pt[0]);
  glEnd();
}

void DrawCircle(Vec2i const&pt, int radius) {
  const int n = 10 + radius / 10;
  glBegin(GL_LINE_STRIP);
  for(int i=0;i<n;++i) {
    Vec2f pi(pt.x, pt.y);
    float a = M_PI * 2 * (float(i) / (n-1));
    pi += Vec2f(cos(a), sin(a)) * radius;
    glVertex2fv(&pi.x);
  }
  glEnd();
}

const int width = 1280, height = 1024;

int random_range(pcg32_random_t *random, int minimum, int maximum) {
  return minimum + pcg32_boundedrand_r(random, maximum - minimum + 1);
}

Vec2i random_edge(pcg32_random_t *random, Extrema2i const&box) {
  switch(random_range(random, 0, 3)) {
    case 0:
      return Vec2i(random_range(random, box.mMin.x, box.mMax.x-1), box.mMin.y);
    case 1:
      return Vec2i(box.mMin.x, random_range(random, box.mMin.y, box.mMax.y-1));
    case 2:
      return Vec2i(random_range(random, box.mMin.x, box.mMax.x-1), box.mMax.y-1);
    case 3:
      return Vec2i(box.mMax.x-1, random_range(random, box.mMin.y, box.mMax.y-1));
  }
  fprintf(stderr, "Shouldn't get here edge");
  exit(1);
}


static std::atomic<unsigned long long> all_steps(0);
static std::atomic<unsigned long long> part_steps_total(0);
static std::atomic<unsigned long long> parts_placed(0);
static std::atomic<int> max_particles(0);
int steps_per_frame = 100;

unique_ptr<pcg32_random_t> random_for_seed(unsigned seed) {
#if 0
  array<uint64_t, 16> seed_ar
          = {{0x374be26ee31f1e78, 0xd4eef394f72f149b, 0x91cb5a7001068c8b, 0x718ef6c2be5efbe7,
              0xbb0dd94396008d70, 0x4f0996d1cd72d2d8, 0x2419b74e0b39e9b3, 0x0da693cf50e1396e,
              0xcaec0e7f4cae7ffa, 0x350b63e4717957c6, 0xbe8460185de680dc, 0xff18c7a0efbcec26,
              0xff1a72bb0ca9ac7f, 0x3b4818e046188158, 0xcac3e320230a44ba, 0xcaf9544740fbd288}};

  for(int i=0;i<16;++i)
    seed_ar[i] += seed;

  unique_ptr<PRNG> ret;
  ret.reset(new PRNG(seed_ar));
#endif
  unique_ptr<pcg32_random_t> ret;
  ret.reset(new pcg32_random_t);
  pcg32_srandom_r(ret.get(), seed, seed >> 16);
  return ret;
}

vector<float> heat;
vector<int> age;

vector<float> burned;

unsigned seed;
Extrema2i aabb;
vector<Vec2i> particles_in_transit;
const int support = 30;
bool run = true;

void trace_heat(Vec2i at) {
  for(size_t iters = 0;(age[at.x + at.y * width] > 0) && (iters < (width*height*10));++iters) {
    assert(age[at.x + at.y * width] >= 0);
    ++heat[at.x + at.y * width];

    int min_neigh_age = age[at.x + at.y * width];

    for(int yy=-1;yy<=1;++yy) {
      for(int xx=-1;xx<=1;++xx) {
        if((xx == 0) && (yy == 0))
          continue;
        const Vec2i neigh = Vec2i(xx, yy) + at;
        if((neigh.x < 0) || (neigh.y < 0))
          continue;
        if((neigh.x >= width) || (neigh.y >= height))
          continue;
        const int neigh_age = age[neigh.x + neigh.y * width];
        if(neigh_age < 0)
          continue;
        if(neigh_age < min_neigh_age) {
          min_neigh_age = neigh_age;
          at = neigh;
        }
      }
    }
  }
}

void bloom_heat() {
  for(int row=0;row<height;++row) {
    for(int col=0;col<width;++col) {

      float total_heat = heat[col + row * width];
      int n_neigh = 1;

      for(int yy=-1;yy<=1;++yy) {
        for(int xx=-1;xx<=1;++xx) {
          if((xx == 0) && (yy == 0))
            continue;
          const Vec2i neigh = Vec2i(xx, yy) + Vec2i(col, row);
          if((neigh.x < 0) || (neigh.y < 0))
            continue;
          if((neigh.x >= width) || (neigh.y >= height))
            continue;
          ++n_neigh;
          const float neigh_heat = heat[neigh.x + neigh.y * width];
          total_heat += neigh_heat / 1.5;
        }
        if(n_neigh > 0)
          heat[col + row * width] = total_heat / n_neigh;
      }
      // Decay
#if 1
      float &heat_here = heat[col + row * width];
      heat_here /= 1.01f;
#endif
    }
  }

  for(int row=0;row<height;++row) {
    for(int col=0;col<width;++col) {
      burned[col + row*width] = std::max(burned[col + row*width], std::min(1.0f, heat[col + row*width]));
    }
  }
}

void OnDraw(void){
	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, width, height, 0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  unique_ptr<pcg32_random_t> random = random_for_seed(rand());

  if(run) {

fprintf(stderr, "go part placements %llu total part-steps %llu total steps %llu steps per placement %f max part %i\n", 
  parts_placed.load(),
  part_steps_total.load(),
  all_steps.load(),
  float(part_steps_total) / parts_placed,
  (int)max_particles.load());

    for(int steps = 0;steps < steps_per_frame;++steps,++all_steps) {
      for(size_t i=0;i<particles_in_transit.size();++i) {
        Vec2i &at = particles_in_transit[i];
        const Vec2i next = 
          at + Vec2i(pcg32_boundedrand_r(random.get(), 3), pcg32_boundedrand_r(random.get(), 3)) - Vec2i(1, 1);

        if((next.x >= 0) && (next.x < width) && (next.y >= 0) && (next.y < height) && 
          (age[next.x + next.y * width] >= 0)) {
          // Fatter
//          if(0 == pcg32_boundedrand_r(random.get(), 16)) 
          {
            aabb.DoEnclose(at - Vec2i(support, support));
            aabb.DoEnclose(at + Vec2i(support, support));
            age[at.x + at.y * width] = age[next.x + next.y * width] + 1;
            trace_heat(at);
            ++parts_placed;

            at = random_edge(random.get(), aabb);
          }
        } else {
          at = next;
          at.x = std::max(aabb.mMin.x, std::min(aabb.mMax.x-1, at.x));
          at.y = std::max(aabb.mMin.y, std::min(aabb.mMax.y-1, at.y));
          ++part_steps_total;
        }
      }

      const int n_particles = ((float(aabb.GetSize().x)+aabb.GetSize().y) * 2000.0f);
      if(particles_in_transit.size() < n_particles) {
        max_particles.store(std::max(max_particles.load(), (int)n_particles));
        particles_in_transit.push_back(random_edge(random.get(), aabb));
      }
      if(0 == (steps % 100))
        bloom_heat();
    }
  }

  vector<float> to_show = burned;

  float min_val = FLT_MAX, max_val = FLT_MIN;
  for(int row=0;row<height;++row) {
    for(int col=0;col<width;++col) {
      min_val = std::min(min_val, to_show[col + row*width]);
      max_val = std::max(max_val, to_show[col + row*width]);
    }
  }

  fprintf(stderr, "display min/max: %f -> %f\n", min_val, max_val);

  vector<unsigned char> buffer;
  buffer.resize(to_show.size());
  for(int row=0;row<height;++row) {
    for(int col=0;col<width;++col) {
      const float norm = ((to_show[col + row*width] - min_val) / (max_val - min_val));
      buffer[col + ((height-1)-row)*width] = 255.0f * (1 - norm);

    }
  }
  for(const auto&pt : particles_in_transit) {
    if((pt.x >= 0) && (pt.x < width) && (pt.y >= 0) && (pt.y < height)) {
      unsigned char &b = buffer[pt.x + ((height-1)-pt.y)*width];
      b = std::min(int(180), int(b) + 64);
    }
  }
  glDrawPixels(width, height, GL_LUMINANCE, GL_UNSIGNED_BYTE, &buffer[0]);

  glColor4f(1,0,0,1);
  DrawBox(aabb);

	glutSwapBuffers();

//  usleep(100000L);
}

void mouseFunc(int button, int state, int x, int y) {
  fprintf(stderr, "mouse %i %i %i\n", button, x, y);
  if(!button && state)
    run = !run;
  else if(button && state)
    steps_per_frame = (steps_per_frame == 1) ? 100 : 1;

  glutPostRedisplay();
}


void idle() {
  glutPostRedisplay();
}



int main(int argc, char** argv) {
  //srand(time(0));
  seed = time(0);

  run = false;

  heat.resize(width*height);
  memset(&heat[0], 0, heat.size() * sizeof(heat[0]));
  burned = heat;

  age.resize(heat.size());
  for(size_t i=0;i<age.size();++i)
    age[i] = -1;

  age[(width/2) + (height/2) * width] = 0;

  aabb = Extrema2i(Vec2i((width/2) - support, (height/2) - support), 
                    Vec2i((width/2) + support, (height/2) + support));

#if 0
  for(int row=0;row<height;++row) {
    for(int col=0;col<width;++col) {
      const int x = col - width/2;
      const int y = row - height/2;
      if(((x*x)+(y*y)) < 100) {
        mheat[col + row*width] = 255;
      }
    }
  }
#endif

  glutInit(&argc, argv);
  glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutCreateWindow("Lichtenberg");

  glutIdleFunc(idle);
	glutDisplayFunc(OnDraw);
  glutMouseFunc(mouseFunc);

	glutMainLoop();
	return 0;
}
