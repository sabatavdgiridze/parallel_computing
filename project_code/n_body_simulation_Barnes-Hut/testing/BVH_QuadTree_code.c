#include "stdlib.h"
#include "math.h"


#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define TRUE 1
#define FALSE 0

struct bounding_box_str {
  double x_min, y_min;
  double x_max, y_max;
};

struct bounding_box_str bounding_box_union(struct bounding_box_str first, struct bounding_box_str second) {
  bounding_box_str union_box;

  union_box.x_min = MIN(first.x_min, second.x_min);
  union_box.x_max = MAX(first.x_max, second.x_max);

  union_box.y_min = MIN(first.y_min, second.y_min);
  union_box.y_max = MAX(first.y_max, second.y_max);

  return union_box;
}


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

struct QUAD_Tree_node {
  struct QUAD_Tree_node * children[4];

  struct particle_str * particle;
  struct bounding_box_str bounding_box;

  double x_center_of_mass;
  double y_center_of_mass;

  double total_mass;
};

struct QUAD_Tree_node_data {
  struct particle_str particle;
  struct bounding_box_str bounding_box;

  double x_center_of_mass;
  double y_center_of_mass;

  double total_mass;
};



struct BVH_node_data {
  struct particle_str particle;
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


void init_bounding_box(struct bounding_box_str * bounding_box, struct particle_str * particle) {
  bounding_box->x_min = particle->x_pos;
  bounding_box->x_max = particle->x_pos;

  bounding_box->y_min = particle->y_pos;
  bounding_box->y_max = particle->y_pos;
}

struct particles_str {
  int n_of_partices;
  struct particle_str * data;
};

double size_bounding_box(struct bounding_box_str bounding_box) {
  double x_span = bounding_box.x_max - bounding_box.x_min;
  double y_span = bounding_box.y_max - bounding_box.y_min;
  return MAX(x_span, y_span);
}

double distance(struct particle_str * particle, bounding_box_str bounding_box) {
  double box_center_x = (bounding_box.x_max + bounding_box.x_min) / 2.0;
  double box_center_y = (bounding_box.y_max + bounding_box.y_min) / 2.0;

  double x_span = box_center_x - particle->x_pos;
  double y_span = box_center_y - particle->y_pos;

  double distance = sqrt(x_span * x_span + y_span * y_span);
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
      max_position = MAX(max_position, current_particle->y_pos);
    } else {
      min_position = MIN(min_position, current_particle->y_pos);
      max_position = MAX(max_position, current_particle->y_pos);
    }
  }
  return max_position - min_position;
}

void swap(struct particle_str * first_particle, struct particle_str * second_particle) {
  struct particle_str * temp = first_particle;
  first_particle = second_particle;
  second_particle = temp;
}

int partition(struct particle_str * particles, int length, int pivot_position, int axis) {
  double pivot_value;
  if (axis == 0) {
    pivot_value = (particles + pivot_position)->x_pos;
  } else {
    pivot_value = (particles + pivot_position)->y_pos;
  }
  swap(particles + pivot_position, particles + length - 1);

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

void select(struct particle_str * particles, int length, int pivot_position, int axis) {
  if (length > 1) {
    int partition_position = rand() % length;
    int pivot_idx = partition(particles, length, partition_position, axis);
    
    if (pivot_position < pivot_idx) {
      select(particles, pivot_idx, pivot_position, axis);
    } else {
      if (pivot_position > pivot_idx) {
        select(particles + (pivot_idx + 1), length - (pivot_idx + 1), pivot_position - (pivot_idx + 1), axis);
      }
    }

  }
}

void swap_double(double * first, double * second) {
  double temporary = *first;
  
  *first = *second;
  *second = temporary;
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

int partition_around_value(struct particle_str * particles, int length, double value, int axis) {
  int left_index = 0;
  int right_index = length ;

  if (axis == 0) {
    while(left_index < right_index) {
      while(left_index < right_index && particles[left_index].x_pos < value) {
        left_index++;
      }
      while(right_index > left_index && particles[right_index - 1].x_pos >= value) {
        right_index--;
      }
      if (left_index < right_index && particles[left_index].x_pos > particles[right_index - 1].x_pos) {
        swap_in_array(particles, left_index, right_index);
        
        left_index++;
        right_index--;
      }
    }
  } else {
    while(left_index < right_index) {
      while(left_index < right_index && particles[left_index].y_pos < value) {
        left_index++;
      }
      while(right_index > left_index && particles[right_index - 1].y_pos >= value) {
        right_index--;
      }
      if (left_index < right_index && particles[left_index].y_pos > particles[right_index - 1].y_pos) {
        swap_in_array(particles, left_index, right_index);

        left_index++;
        right_index--;
      }
    }
  }
  return right_index;
}

struct bounding_box_str * divide_bounding_box(struct bounding_box_str bounding_box) {
  struct bounding_box_str * boxes = (bounding_box_str *) malloc(sizeof(struct bounding_box_str) * 4);
  
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



struct QUAD_Tree_node * create_QUAD_tree(struct particle_str * particles, struct bounding_box_str bounding_box, int length) {
  if (length == 0 || length == 1) {
    struct QUAD_Tree_node * node = (struct QUAD_Tree_node *) malloc(sizeof(struct QUAD_Tree_node));

    int child_idx;
    for (child_idx = 0; child_idx < 4; child_idx++) {
      node->children[child_idx] = NULL;
    }
    node->particle = particles;
    node->bounding_box = bounding_box;
    if (particles != NULL) {
      node->x_center_of_mass = particles->x_pos;
      node->y_center_of_mass = particles->y_pos;

      node->total_mass = particles->mass;
    }
    else {
      node->x_center_of_mass = 0.0;
      node->y_center_of_mass = 0.0;

      node->total_mass = 0.0;
    }
    return node;
  } else {
    struct QUAD_Tree_node * node = (struct QUAD_Tree_node *) malloc(sizeof(struct QUAD_Tree_node));

    node->particle = NULL;
    node->bounding_box = bounding_box;
    
    struct bounding_box_str * boxes = divide_bounding_box(bounding_box);


    double horizontal_division = (bounding_box.x_max + bounding_box.x_min) / 2.0;
    double vertical_division = (bounding_box.y_max + bounding_box.y_min) / 2.0;

    int horizontal_pivot = partition_around_value(particles, length, horizontal_division, 0);
    if (horizontal_pivot > 0) {
      // left side
      
      int left_length = horizontal_pivot;
      struct particle_str * left_particles = particles;


      int vertical_pivot = partition_around_value(left_particles, left_length, vertical_division, 1);

      if (vertical_pivot > 0) {
        node->children[1] = create_QUAD_tree(left_particles, boxes[1], vertical_pivot);
      } else {
        node->children[1] = create_QUAD_tree(NULL, boxes[1], 0);
      }
      if (vertical_pivot < left_length) {
        node->children[0] = create_QUAD_tree(left_particles + vertical_pivot, boxes[0], left_length - vertical_pivot);
      } else {
        node->children[0] = create_QUAD_tree(NULL, boxes[0], 0);
      }

    } else {
      // left side is empty
      node->children[0] = create_QUAD_tree(NULL, boxes[0], 0);
      node->children[1] = create_QUAD_tree(NULL, boxes[1], 0);
    }
    if (horizontal_pivot < length) {
      // right side

      int right_length = length - horizontal_pivot;
      struct particle_str * right_particles = particles + horizontal_pivot;


      int vertical_pivot = partition_around_value(right_particles, right_length, vertical_division, 1);

      if (vertical_pivot > 0) {
        node->children[2] = create_QUAD_tree(right_particles, boxes[2], vertical_pivot);
      } else {
        node->children[2] = create_QUAD_tree(NULL, boxes[2], 0);
      }
      if (vertical_pivot < right_length) {
        node->children[3] = create_QUAD_tree(right_particles + vertical_pivot, boxes[3], right_length - vertical_pivot);
      } else {
        node->children[3] = create_QUAD_tree(NULL, boxes[3], 0);
      }

    } else {
      // right side is empty
      node->children[2] = create_QUAD_tree(NULL, boxes[2], 0);
      node->children[3] = create_QUAD_tree(NULL, boxes[3], 0);
    }
    int child_idx;
    
    node->total_mass = 0.0;
    node->x_center_of_mass = 0.0;
    node->y_center_of_mass = 0.0;

    for (child_idx = 0; child_idx < 4; child_idx++) {
      node->total_mass += node->children[child_idx]->total_mass;
      
      node->x_center_of_mass += node->children[child_idx]->x_center_of_mass * node->children[child_idx]->total_mass;
      node->y_center_of_mass += node->children[child_idx]->y_center_of_mass * node->children[child_idx]->total_mass;
    }
    node->x_center_of_mass /= node->total_mass;
    node->y_center_of_mass /= node->total_mass;

    return node;
  }
}

struct BVH_node * create_BVH(struct particle_str * particles,int length) {
  if (length == 1) {
    BVH_node * node = (BVH_node *) malloc(sizeof(struct BVH_node));
    
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
    select(particles, length, first_child_length - 1, division_axis);

    BVH_node * first_child = create_BVH(particles, first_child_length);
    BVH_node * second_child = create_BVH(particles + first_child_length, length - first_child_length);


    BVH_node * node = (BVH_node *) malloc(sizeof(struct BVH_node));



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


void flatten_QUAD_Tree(struct QUAD_Tree_node * root, struct QUAD_Tree_node_data * particles_data, int * first_child, int * next_sibling) {
  // initially, both first_child and next_sibling should be assigned -1
  
  flatten_QUAD_Tree_rec(0, root, particles_data, first_child, next_sibling);
}


int flatten_QUAD_Tree_rec(int position, struct QUAD_Tree_node * root, struct QUAD_Tree_node_data * particles_data, int * first_child, int * next_sibling) {
  (particles_data + position)->bounding_box = root->bounding_box;
  (particles_data + position)->particle = *(root->particle);

  (particles_data + position)->x_center_of_mass = root->x_center_of_mass;
  (particles_data + position)->y_center_of_mass = root->y_center_of_mass;
  (particles_data + position)->total_mass = root->total_mass;

  if (root->children[0] != NULL) {
    first_child[position] = position + 1;
  } 
  
  int child_idx;
  int current_position = position + 1;
  for (child_idx = 0; child_idx < 4; child_idx++) {
    int next_position = flatten_QUAD_Tree_rec(current_position, root->children[child_idx], particles_data, first_child, next_sibling);
    if (child_idx < 4) {
      next_sibling[current_position] = next_position;
    }
    current_position = next_position;
  }
}

void flatten_BVH_tree(struct BVH_node * root, struct BVH_node_data * particles_data, int * next) {
  flatten_tree_BVH_rec(0, root, particles_data, next);
}


int flatten_tree_BVH_rec(int position, struct BVH_node * root, struct BVH_node_data * particles_data, int * next) {

  (particles_data + position)->bounding_box = root->bounding_box;
  (particles_data + position)->particle = *(root->particle);

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


void compute_force(struct particle_str * particle, struct particle_str * other, double gravity_constant) {
  double x_sep, y_sep;
  
  x_sep = other->x_pos - particle->x_pos;
  y_sep = other->y_pos - particle->y_pos;

  double distance_squared = x_sep * x_sep + y_sep * y_sep;

  double grav_base = gravity_constant * (particle->mass) * (other->mass) / distance_squared;

  particle->x_force += grav_base * x_sep;
  particle->y_force += grav_base * y_sep;
}


void compute_force_BVH_tree(struct particle_str * particle, BVH_node * root, double gravity_constant, double threshold) {
  if (root->particle == NULL) {
    if (size_bounding_box(root->bounding_box) / distance(particle, root->bounding_box) < threshold) {
      struct particle_str * particle_union = (particle_str *) malloc(sizeof(struct particle_str));
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

void compute_force_QUAD_Tree(struct particle_str * particle, struct QUAD_Tree_node * root, double gravity_constant, double threshold) {
  if (root->particle == NULL) {
    if (root->children[0] != NULL) {
      // node is non-empty
      if (size_bounding_box(root->bounding_box) / distance(particle, root->bounding_box) < threshold) {
        struct particle_str * particle_union = (struct particle_str *) malloc(sizeof(struct particle_str));
        particle_union->mass = root->total_mass;

        particle_union->x_pos = root->x_center_of_mass;
        particle_union->y_pos = root->y_center_of_mass;

        compute_force(particle, particle_union, gravity_constant);
      } else {
        int child_idx;
        for (child_idx = 0; child_idx < 4; child_idx++) {
          compute_force_QUAD_Tree(particle, root->children[child_idx], gravity_constant, threshold);
        }
      }
    }
  } else {
    compute_force(particle, root->particle, gravity_constant);
  }
}


void compute_force_flattened_BVH_tree(struct particle_str * particle, struct BVH_node_data * particles_data, int * next, double gravity_constant, double threshold) {
  int max_depth;
  int * node_to_visit = (int *) malloc(sizeof(int) * max_depth);

  int to_visit_idx = 0;
  int current_node_idx = 0;

  while (TRUE) {
    if (size_bounding_box((particles_data + current_node_idx)->bounding_box) / distance(particle, (particles_data + current_node_idx)->bounding_box) < threshold) {
      struct particle_str * particle_union = (particle_str *) malloc(sizeof(struct particle_str));
      
      particle_union->mass = (particles_data + current_node_idx)->total_mass;

      particle_union->x_pos = (particles_data + current_node_idx)->x_center_of_mass;
      particle_union->y_pos = (particles_data + current_node_idx)->y_center_of_mass;

      compute_force(particle, particle_union, gravity_constant);

      if (to_visit_idx == 0) {
        break;
      } else {
        current_node_idx = node_to_visit[--to_visit_idx];
      }
    } else {
      node_to_visit[to_visit_idx++] = current_node_idx + 1;
      node_to_visit[to_visit_idx++] = next[current_node_idx];
    }
  }
}