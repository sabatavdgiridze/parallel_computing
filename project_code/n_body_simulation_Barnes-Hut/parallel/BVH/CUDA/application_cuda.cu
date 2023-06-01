#include <cuda_runtime.h>
#include <cuda.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define TRUE 1
#define FALSE 0

#define LOG_THE_DATA

//#define DEBUG

#define NO_DEBUG


struct bounding_box_str {
  double x_min, y_min;
  double x_max, y_max;
};

struct BVH_node {
  struct BVH_node * first_child;
  struct BVH_node * second_child;

  int division_axis;

  struct particle_str  * particle;
  struct bounding_box_str bounding_box;

  double x_center_of_mass;
  double y_center_of_mass;

  double total_mass;
};


struct particle_str {
  double x_pos, y_pos;
  double x_vel, y_vel;
  double x_force, y_force;
  double mass;
};


struct BVH_node_data {
  struct bounding_box_str bounding_box;

  double x_center_of_mass;
  double y_center_of_mass;

  double total_mass;
};



struct bounding_box_str bounding_box_union(struct bounding_box_str first, struct bounding_box_str second) {
  struct bounding_box_str union_box;

  union_box.x_min = MIN(first.x_min, second.x_min);
  union_box.x_max = MAX(first.x_max, second.x_max);

  union_box.y_min = MIN(first.y_min, second.y_min);
  union_box.y_max = MAX(first.y_max, second.y_max);

  return union_box;
}


void init_bounding_box(struct bounding_box_str * bounding_box, struct particle_str * particle) {
  bounding_box->x_min = particle->x_pos;
  bounding_box->x_max = particle->x_pos;

  bounding_box->y_min = particle->y_pos;
  bounding_box->y_max = particle->y_pos;
}

double size_bounding_box(struct bounding_box_str bounding_box) {
  double x_span = bounding_box.x_max - bounding_box.x_min;
  double y_span = bounding_box.y_max - bounding_box.y_min;
  return MAX(x_span, y_span);
}

double distance(struct particle_str * particle, struct bounding_box_str bounding_box) {
  double box_center_x = (bounding_box.x_max + bounding_box.x_min) / 2.0;
  double box_center_y = (bounding_box.y_max + bounding_box.y_min) / 2.0;

  double x_span = box_center_x - particle->x_pos;
  double y_span = box_center_y - particle->y_pos;

  double distance = sqrt(x_span * x_span + y_span * y_span);
  return distance;
}

double calculate_variance(struct particle_str * particles, int length, int axis) {
  double min_position;
  double max_position;
  
  if (axis == 0) {
    min_position = particles->x_pos;
    max_position = particles->x_pos;
  } else {
    min_position = particles->y_pos;
    max_position = particles->y_pos;
  }
  int idx;
  for (idx = 1; idx < length; idx++) {
    struct particle_str * current_particle = particles + idx;
    if (axis == 0) {
      min_position = MIN(min_position, current_particle->x_pos);
      max_position = MAX(max_position, current_particle->x_pos);
    } else {
      min_position = MIN(min_position, current_particle->y_pos);
      max_position = MAX(max_position, current_particle->y_pos);
    }
  }
  return max_position - min_position;
}


void swap_double(double * first, double * second) {
    double temporary = * first;
    * first = * second;
    * second = temporary;
}




void swap(struct particle_str * first_particle, struct particle_str * second_particle) {
  swap_double(&(first_particle->x_pos), &(second_particle->x_pos));
  swap_double(&(first_particle->y_pos), &(second_particle->y_pos));

  swap_double(&(first_particle->x_vel), &(second_particle->x_vel));
  swap_double(&(first_particle->y_vel), &(second_particle->y_vel));
 
  swap_double(&(first_particle->x_force), &(second_particle->x_force));
  swap_double(&(first_particle->y_force), &(second_particle->y_force));

  swap_double(&(first_particle->mass), &(second_particle->mass));
}

int partition(struct particle_str * particles, int length, int pivot_position, int axis) {
  double pivot_value;
  if (axis == 0) {
    pivot_value = (particles + pivot_position)->x_pos;
  } else {
    pivot_value = (particles + pivot_position)->y_pos;
  }
  swap(&particles[pivot_position], &particles[length - 1]);


  int pivot_idx = 0;
  int idx;
  for (idx = 0; idx < length - 1; idx++) {
    if (axis == 0) {
      if ((particles+idx)->x_pos < pivot_value) {
        swap(particles + pivot_idx, particles + idx);
        pivot_idx++;
      }
    } else {
      if ((particles+idx)->y_pos < pivot_value) {
        swap(particles + pivot_idx, particles + idx);
        pivot_idx++;
      }
    }
  }
  swap(particles + pivot_idx, particles + length - 1);
  return pivot_idx;
}


void select_BVH(struct particle_str * particles, int length, int pivot_position, int axis) {
  if (length > 1) {
    int partition_position = rand() % length;
    int pivot_idx = partition(particles, length, partition_position, axis);
    
    if (pivot_position < pivot_idx) {
      select_BVH(particles, pivot_idx, pivot_position, axis);
    } else {
      if (pivot_position > pivot_idx) {
        select_BVH(particles + (pivot_idx + 1), length - (pivot_idx + 1), pivot_position - (pivot_idx + 1), axis);
      }
    }

  }
}

void swap_in_array(struct particle_str * particles, int left_index, int right_index) {
  swap_double(&particles[left_index].mass, &particles[right_index].mass);

  swap_double(&particles[left_index].x_pos, &particles[right_index].x_pos);
  swap_double(&particles[left_index].y_pos, &particles[right_index].y_pos);
  
  swap_double(&particles[left_index].x_vel, &particles[right_index].x_vel);
  swap_double(&particles[left_index].y_vel, &particles[right_index].y_vel);

  swap_double(&particles[left_index].x_force, &particles[right_index].x_force);
  swap_double(&particles[left_index].y_force, &particles[right_index].y_force);
}

struct bounding_box_str * divide_bounding_box(struct bounding_box_str bounding_box) {
  struct bounding_box_str * boxes = (struct bounding_box_str *) malloc(sizeof(struct bounding_box_str) * 4);
  
  double x_span = (bounding_box.x_max - bounding_box.x_min) / 2.0;
  double y_span = (bounding_box.y_max - bounding_box.y_min) / 2.0;

  double x_offset[4] = {0, 0, 1, 1};
  double y_offset[4] = {1, 0, 0, 1};

  int box_idx;
  for (box_idx = 0; box_idx < 4; box_idx++) {
    boxes[box_idx].x_min = bounding_box.x_min + x_offset[box_idx] * x_span;
    boxes[box_idx].y_min = bounding_box.y_min + y_offset[box_idx] * y_span;

    boxes[box_idx].x_max = boxes[box_idx].x_min + x_span;
    boxes[box_idx].y_max = boxes[box_idx].y_min + y_span;
  }

  return boxes;
}

struct BVH_node * memory;
int stack_pointer;


void init_memory() {
    int N = 10000000;
    stack_pointer = 0;
    memory = (struct BVH_node *) malloc(sizeof(struct BVH_node) * N);
}

struct BVH_node * create_node() {
    int old_stack_pointer = stack_pointer;
    stack_pointer++;

    return (memory + old_stack_pointer);
}




struct BVH_node * create_BVH(struct particle_str * particles,int length) {  
  if (length == 1) {
    struct BVH_node * node = create_node();
    
    node->first_child = NULL;
    node->second_child = NULL;
    
    node->division_axis = -1;
    node->particle = particles;

    init_bounding_box(&(node->bounding_box), particles);

    node->x_center_of_mass = particles->x_pos;
    node->y_center_of_mass = particles->y_pos;

    node->total_mass = particles->mass;

    return node;
  } else {
    int division_axis;
    
    double variance_x = calculate_variance(particles, length, 0);
    double variance_y = calculate_variance(particles, length, 1);

    

    if (variance_x > variance_y) {
      division_axis = 0;
    } else {
      division_axis = 1;
    }

    int first_child_length = length / 2;
    select_BVH(particles, length, first_child_length - 1, division_axis);

    struct BVH_node * first_child = create_BVH(particles, first_child_length);
    struct BVH_node * second_child = create_BVH(particles + first_child_length, length - first_child_length);

    struct BVH_node * node = create_node();


    node->first_child = first_child;
    node->second_child = second_child;

    node->division_axis = division_axis;
    node->particle = NULL;

    node->bounding_box = bounding_box_union(first_child->bounding_box, second_child->bounding_box);

    node->total_mass = first_child->total_mass + second_child->total_mass;
    node->x_center_of_mass = (first_child->x_center_of_mass * first_child->total_mass + second_child->x_center_of_mass * second_child->total_mass) / node->total_mass;
    node->y_center_of_mass = (first_child->y_center_of_mass * first_child->total_mass + second_child->y_center_of_mass * second_child->total_mass) / node->total_mass;
    
    return node;
  }
}

void compute_force(struct particle_str * particle, struct particle_str * other, double gravity_constant) {
  double x_sep, y_sep;
  
  x_sep = other->x_pos - particle->x_pos;
  y_sep = other->y_pos - particle->y_pos;

  double distance_squared = x_sep * x_sep + y_sep * y_sep;

  double grav_base = gravity_constant * (particle->mass) * (other->mass) / distance_squared;

  particle->x_force += grav_base * x_sep;
  particle->y_force += grav_base * y_sep;
}


void compute_force_BVH_tree(struct particle_str * particle, struct BVH_node * root, double gravity_constant, double threshold) {
  if (root->particle == NULL) {
    if (size_bounding_box(root->bounding_box) / distance(particle, root->bounding_box) < threshold) {
      struct particle_str * particle_union = (struct particle_str *) malloc(sizeof(struct particle_str));
      particle_union->mass = root->total_mass;

      particle_union->x_pos = root->x_center_of_mass;
      particle_union->y_pos = root->y_center_of_mass;

      compute_force(particle, particle_union, gravity_constant);

    } else {
      compute_force_BVH_tree(particle, root->first_child, gravity_constant, threshold);
      compute_force_BVH_tree(particle, root->second_child, gravity_constant, threshold);
    }
  } else {
    compute_force(particle, root->particle, gravity_constant);
  }
}


void print_indentation(int indentation) {
    int idx;
    for(idx = 0; idx < indentation; idx++) {
        printf(" ");
    }
}


void pretty_print_BVH(struct BVH_node * root) {
    static int indentation = 0;
    if (root->particle != NULL) {
        print_indentation(indentation);
        printf("particle = (%f, %f)\n", root->particle->x_pos, root->particle->y_pos);
    } else {

        print_indentation(indentation);
        if (root->division_axis == 0) {
            printf("divided on x-axis\n");
        } else {
            printf("divided on y-axis\n");
        }
        print_indentation(indentation);
        printf("bound = (x-min = %f, x-max = %f, y-min = %f, y-max = %f)\n\n", root->bounding_box.x_min, root->bounding_box.x_max, root->bounding_box.y_min, root->bounding_box.y_max);

        indentation += 4;
        pretty_print_BVH(root->first_child);
        
        pretty_print_BVH(root->second_child);

        indentation -= 4;
    }
}

int flatten_tree_BVH_rec(int position, struct BVH_node * root, struct BVH_node_data * particles_data, int * next);

void flatten_BVH_tree(struct BVH_node * root, struct BVH_node_data * particles_data, int * next) {
    flatten_tree_BVH_rec(0, root, particles_data, next);
}


int flatten_tree_BVH_rec(int position, struct BVH_node * root, struct BVH_node_data * particles_data, int * next) {

    if (root->particle == NULL) {
      (particles_data + position)->bounding_box = root->bounding_box;
    } else {
      init_bounding_box(&((particles_data + position)->bounding_box), root->particle);
    }
    (particles_data + position)->x_center_of_mass = root->x_center_of_mass;
    (particles_data + position)->y_center_of_mass = root->y_center_of_mass;
    (particles_data + position)->total_mass = root->total_mass;
    
    if (root->first_child == NULL) {
        next[position] = -1;
        return position + 1;
    } else {
        next[position] = flatten_tree_BVH_rec(position + 1, root->first_child, particles_data, next);
        return flatten_tree_BVH_rec(next[position], root->second_child, particles_data, next);
    }
}

void compute_force_flattened_BVH_tree(struct particle_str * particle, struct BVH_node_data * particles_data, int * next, double gravity_constant, double threshold) {
  int max_depth = 1000;
  int * node_to_visit = (int *) malloc(sizeof(int) * max_depth);

  int to_visit_idx = 0;
  int current_node_idx = 0;
  node_to_visit[to_visit_idx] = 0;

  while (TRUE) {
    current_node_idx = node_to_visit[to_visit_idx];

    if (size_bounding_box((particles_data + current_node_idx)->bounding_box) / distance(particle, (particles_data + current_node_idx)->bounding_box) < threshold) {
      struct particle_str * particle_union = (struct particle_str *) malloc(sizeof(struct particle_str));
      
      particle_union->mass = (particles_data + current_node_idx)->total_mass;

      particle_union->x_pos = (particles_data + current_node_idx)->x_center_of_mass;
      particle_union->y_pos = (particles_data + current_node_idx)->y_center_of_mass;

      compute_force(particle, particle_union, gravity_constant);

      if (to_visit_idx == 0) {
        break;
      } else {
        to_visit_idx--;
      }
    } else {
      if (next[current_node_idx] != -1) {
        node_to_visit[to_visit_idx++] = current_node_idx + 1;
        node_to_visit[to_visit_idx] = next[current_node_idx];
      }
    }
  }
}



#define MAX_DEPTH 1000

__global__ void cuda_flatten_tree_update(int N_array, double * device_particles_array, double * device_flattened_tree_array, int * device_next, double gravity_constant, double threshold, double time_step, double mass) {
    int idx = blockDim.x * blockIdx.x + threadIdx.x;

    int node_to_visit[MAX_DEPTH];
    while (idx < N_array) {
        int to_visit_idx = 0;
        int current_node_idx = 0;
        node_to_visit[to_visit_idx] = current_node_idx;

        while (TRUE) {
            current_node_idx = node_to_visit[to_visit_idx];
            double size_bounding_box_x = device_flattened_tree_array[7 * current_node_idx + 1] - device_flattened_tree_array[7 * current_node_idx    ];
            double size_bounding_box_y = device_flattened_tree_array[7 * current_node_idx + 3] - device_flattened_tree_array[7 * current_node_idx + 2];

            double size_bounding_box = MAX(size_bounding_box_x, size_bounding_box_y);

            double center_x = device_flattened_tree_array[7 * current_node_idx + 4];
            double center_y = device_flattened_tree_array[7 * current_node_idx + 5];

            double pos_x = device_particles_array[6 * idx    ];
            double pos_y = device_particles_array[6 * idx + 1];

            double distance_squared = (center_x - pos_x) * (center_x - pos_x) + (center_y - pos_y) * (center_y - pos_y);
            double eps = 0.0001;

            if ( (size_bounding_box * size_bounding_box) < distance_squared  < threshold) {
                double union_mass = device_flattened_tree_array[7 * current_node_idx + 6];

                double grav_base = gravity_constant * mass * union_mass / distance_squared;

                device_particles_array[6 * idx + 4] += grav_base * (center_x - pos_x);
                device_particles_array[6 * idx + 5] += grav_base * (center_y - pos_y);

                if (to_visit_idx == 0) {
                    break;
                } else {
                    to_visit_idx--;
                }
            } else {
                if (device_next[current_node_idx] != -1) {
                    node_to_visit[to_visit_idx++] = current_node_idx + 1;
                    node_to_visit[to_visit_idx] = device_next[current_node_idx];
                }
            }
        }

        idx += gridDim.x * blockDim.x;
    }
    
    idx = blockDim.x * blockIdx.x + threadIdx.x;
    while(idx < N_array) {
        device_particles_array[6 * idx + 2] += device_particles_array[6 * idx + 4] / mass * time_step;
        device_particles_array[6 * idx + 3] += device_particles_array[6 * idx + 5] / mass * time_step;

        device_particles_array[6 * idx    ] += device_particles_array[6 * idx + 2] * time_step;
        device_particles_array[6 * idx + 1] += device_particles_array[6 * idx + 3] * time_step;
        
        
        idx += gridDim.x * blockDim.x;
    }
}
void compute_force_and_update_particles_flattened_BVH_tree_CUDA(int N_array, struct particle_str * particles, struct BVH_node_data * flattened_BVH_tree, int * next, double gravity_constant, double threshold, double time_step) {
    double * particles_array = (double *) malloc(sizeof(double) *  6 * N_array);
    int idx ;

    double mass = particles[0].mass;

    for (idx = 0; idx < N_array; idx++) {
        particles_array[6 * idx    ] = particles[idx].x_pos;
        particles_array[6 * idx + 1] = particles[idx].y_pos;
        particles_array[6 * idx + 2] = particles[idx].x_vel;
        particles_array[6 * idx + 3] = particles[idx].y_vel;
        particles_array[6 * idx + 4] = particles[idx].x_force;
        particles_array[6 * idx + 5] = particles[idx].y_force;
    }
    double * flattened_tree_array = (double *) malloc(sizeof(double) *  7 * 4 * N_array);
    for (idx = 0; idx < 4 * N_array; idx++) {

        flattened_tree_array[7 * idx    ] = flattened_BVH_tree[idx].bounding_box.x_min;
        flattened_tree_array[7 * idx + 1] = flattened_BVH_tree[idx].bounding_box.x_max;
        flattened_tree_array[7 * idx + 2] = flattened_BVH_tree[idx].bounding_box.y_min;
        flattened_tree_array[7 * idx + 3] = flattened_BVH_tree[idx].bounding_box.y_max;
        flattened_tree_array[7 * idx + 4] = flattened_BVH_tree[idx].x_center_of_mass;
        flattened_tree_array[7 * idx + 5] = flattened_BVH_tree[idx].y_center_of_mass;
        flattened_tree_array[7 * idx + 6] = flattened_BVH_tree[idx].total_mass;
    }

    double * device_particles_array;
    double * device_flattened_tree_array;
    int * device_next;

    cudaMalloc((void**)&device_particles_array, 6 * N_array * sizeof(double));
    cudaMalloc((void**)&device_flattened_tree_array, 7 * 4 *N_array * sizeof(double));
    cudaMalloc((void**)&device_next, 4 * N_array * sizeof(int));

    cudaMemcpy(device_particles_array, particles_array, 6 * N_array * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_flattened_tree_array, flattened_tree_array, 7 * 4 * N_array * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_next, next, 4 * N_array * sizeof(int), cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (N_array + threadsPerBlock - 1) / threadsPerBlock;

    cuda_flatten_tree_update<<<blocksPerGrid, threadsPerBlock>>>(N_array, device_particles_array, device_flattened_tree_array, device_next, gravity_constant, threshold, time_step, mass);


    cudaMemcpy(particles_array, device_particles_array, 6 * N_array * sizeof(double), cudaMemcpyDeviceToHost);
    for (idx = 0; idx < N_array; idx++) {
        particles[idx].x_pos = particles_array[6 * idx    ];
        particles[idx].y_pos = particles_array[6 * idx + 1];
        particles[idx].x_vel = particles_array[6 * idx + 2];
        particles[idx].y_vel = particles_array[6 * idx + 3];
        particles[idx].x_force = 0.0;
        particles[idx].y_force = 0.0;
    }
    free(particles_array);
    free(flattened_tree_array);
    cudaFree(device_particles_array);
    cudaFree(device_flattened_tree_array);
    cudaFree(device_next);
}




double get_random() {
    return (1.0 * rand() / RAND_MAX) * 2 - 1.0;
}

struct particle_str * init_particles(int N_array) {
    struct particle_str * particles = (struct particle_str *) malloc(sizeof(struct particle_str) * N_array);

    double MASS = 0.1;

    int idx;
    for (idx = 0; idx < N_array; idx++) {
        particles[idx].mass = MASS;

        particles[idx].x_pos = get_random();
        particles[idx].y_pos = get_random();

        particles[idx].x_vel = get_random() * 0.1;
        particles[idx].y_vel = get_random() * 0.1;

        particles[idx].x_force = 0.0;
        particles[idx].y_force = 0.0;

    }
    return particles;
}

void update_particles(struct particle_str * particles, int N_array, double delta_t) {
    int idx;
    for (idx = 0; idx < N_array; idx++) {
        particles[idx].x_vel += particles[idx].x_force / particles[idx].mass * delta_t;
        particles[idx].y_vel += particles[idx].y_force / particles[idx].mass * delta_t;

        particles[idx].x_pos += particles[idx].x_vel * delta_t;
        particles[idx].y_pos += particles[idx].y_vel * delta_t;

        particles[idx].x_force = 0.0;
        particles[idx].y_force = 0.0;
    }
}



void application(int N_array) {
    init_memory();
    struct particle_str * particles = init_particles(N_array);
    
    struct BVH_node_data * flattened_BVH_tree = (struct BVH_node_data *) malloc(sizeof(struct BVH_node_data) * 4 * N_array);
    int * next = (int *) malloc(sizeof(int) * 4 * N_array);
    
    int next_idx;
    for (next_idx = 0; next_idx < 4 * N_array; next_idx++) {
      next[next_idx] = -1;
    }


    double time_start = 0.0;
    double time_finish = 1.0;

    double time_step = 0.01;

    double time;

    double gravity_constant = 0.1;
    double threshold = 2.0;

    clock_t timer_start;
    clock_t timer_end;


    timer_start = clock();
    for (time = time_start; time <= time_finish; time++) {
        struct BVH_node * BVH_root = create_BVH(particles, N_array);
        
        flatten_BVH_tree(BVH_root, flattened_BVH_tree, next);

        int idx;
        for (idx = 0; idx < N_array; idx++) {

            compute_force_and_update_particles_flattened_BVH_tree_CUDA(N_array, particles, flattened_BVH_tree, next, gravity_constant, threshold, time_step);

        }
        stack_pointer = 0;
    }
    timer_end = clock();
#ifdef DEBUG
    int print_idx;
    for (print_idx = 0; print_idx < N_array; print_idx++) {
        printf("(%f, %f) ", particles[print_idx].x_pos, particles[print_idx].y_pos);
    }
    printf("\n");
#endif
    printf("It took : %d milliseconds!\n", (timer_end - timer_start) * 1000 / CLOCKS_PER_SEC);
}




int main(int argc, char ** argv) {
    int N_array = atoi(argv[1]);
    int Seed = atoi(argv[2]);

    srand(Seed);
    application(N_array);

    return 0;
}