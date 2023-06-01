#include <cuda_runtime.h>
#include <cuda.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

#define STACK_MAX 100

int color_count = 6;

__global__ void compute_num(int * queries, int queries_length, int * keys, int * has_children, int * next_sibling, int * answer) {
  int stack[STACK_MAX];

  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  while (idx < queries_length) {
    int stack_pointer = 0;
    stack[stack_pointer] = 0;

    while (TRUE) {
      if (keys[stack[stack_pointer]] == queries[idx]) {
        atomicAdd(&(answer[idx]), 1);
      }

      if (has_children[stack[stack_pointer]] == FALSE) {
        stack_pointer--;
      } else {
        int child_pos = stack[stack_pointer] + 1;

        stack_pointer--;
        while(child_pos != -1) {
          stack_pointer++;
          stack[stack_pointer] = child_pos;
          child_pos = next_sibling[child_pos];
        }
      }
      if (stack_pointer == -1) {
        break;
      }
    }

    idx += gridDim.x * blockDim.x;
  }
}


void cuda_compute(int N, int * keys, int * has_children, int * next_position) {
  int * device_keys;
  int * device_has_children;
  int * device_next_position;

  cudaMalloc(&device_keys, sizeof(int) * N);
  cudaMemcpy(device_keys, keys, sizeof(int) * N, cudaMemcpyHostToDevice);

  cudaMalloc(&device_has_children, sizeof(int) * N);
  cudaMemcpy(device_has_children, has_children, sizeof(int) * N, cudaMemcpyHostToDevice);

  cudaMalloc(&device_next_position, sizeof(int) * N);
  cudaMemcpy(device_next_position, next_position, sizeof(int) * N, cudaMemcpyHostToDevice);


  


  int queries_length = color_count;
  int * queries = (int *) malloc(sizeof(int) * queries_length);
  int * answers = (int *) malloc(sizeof(int) * queries_length);


  int idx;
  for (idx = 0; idx < queries_length; idx++) {
    queries[idx] = idx;
  }

  int threadsPerBlock = 1;
  int blocksPerGrid = 1;

  int * device_queries;
  int * device_answers;
  
  cudaMalloc(&device_queries, sizeof(int) * queries_length);
  cudaMemcpy(device_queries, queries, sizeof(int) * queries_length, cudaMemcpyHostToDevice);
  
  cudaMalloc(&device_answers, sizeof(int) * queries_length);

  compute_num<<<blocksPerGrid, threadsPerBlock>>>(device_queries, queries_length, device_keys, device_has_children, device_next_position, device_answers);
  cudaMemcpy(answers, device_answers, sizeof(int) * queries_length, cudaMemcpyDeviceToHost);

  for (idx = 0; idx < queries_length; idx++) {
    printf("(%d -> %d)\n", idx, answers[idx]);
  }

}

struct node_str {
  int key;

  int child_count;
  struct node_str ** children;
};

struct node_str * create_node(int key) {
  struct node_str * node = (struct node_str *) malloc(sizeof(struct node_str));
  node->key = key;

  node->child_count = 0;
  node->children = NULL;

  return node;
}

void add_node(struct node_str * root, int key) {
  root->child_count++;
  root->children = (struct node_str **) realloc(root->children, (root->child_count) * sizeof(struct node_str *));
  root->children[root->child_count - 1] = create_node(key);
}


int get_number(int low, int high) {
  return (low + rand() % (high - low + 1));
}

void print_tree(struct node_str * root);

struct node_str * create_tree(int N) {
  struct node_str * root = create_node(get_number(0, color_count - 1));
  
  int nodes_count = 1;

  struct node_str * stack[STACK_MAX];

  int stack_pointer = 0;
  stack[stack_pointer++] = root;

  while (nodes_count < N) {

    int children_min;
    int children_max;
    if (stack_pointer == 1) {
      children_min = 1;
    } else {
      children_min = 0;
    }


    struct node_str * current_node = stack[--stack_pointer];


    children_max = N - nodes_count;
    int child_count = get_number(children_min, children_max);
    nodes_count += child_count;

    int child_idx;
    for (child_idx = 0; child_idx < child_count; child_idx++) {
      add_node(current_node, get_number(0, color_count - 1));
      stack[stack_pointer++] = current_node->children[child_idx];
    }
  }
  return root;
}


int flatten_tree_rec(int position, struct node_str * root, int * keys, int * has_children, int * next_sibling) {
  keys[position] = root->key;
  if (root->child_count == 0) {
    has_children[position]= FALSE;
    return position + 1;
  }
  else {

    has_children[position] = TRUE;

    position++;
    int next_position;

    struct node_str * current_node;
    int child_idx;
    for (child_idx = 0; child_idx < root->child_count; child_idx++) {
      current_node = root->children[child_idx];
      next_position = flatten_tree_rec(position, current_node, keys, has_children, next_sibling);

      if (child_idx < root->child_count - 1) {
        next_sibling[position] = next_position;
      }
      if (child_idx == root->child_count - 1) {
        next_sibling[position] = -1;
      }
      position = next_position;
    }
    return position;
  }
}


void flatten_tree(struct node_str * root, int * keys, int * has_children, int * next_sibling) {
  flatten_tree_rec(0, root, keys, has_children, next_sibling);
}

void print_character(char c, int times) {
  int idx;
  for (idx = 0; idx < times; idx++) {
    printf("%c", c);
  }
}

void print_tree(struct node_str * root) {
  static int indentation = 0;
  print_character(' ', indentation);
  printf("%d\n", root->key);

  if(root->child_count > 0) {
    int child_idx;
    for (child_idx = 0; child_idx < root->child_count; child_idx++) {
      indentation += 2;
      print_tree(root->children[child_idx]);
      indentation -= 2;
    }
  }
}

void print_flattened_tree(int root, int * keys, int * has_children, int * next_sibling) {
  static int indentation = 0;
  print_character(' ', indentation);
  printf("%d\n", keys[root]);
  if (has_children[root] == TRUE) {
    int next_position = root + 1;
    while(next_position != -1) {
      indentation += 2;
      print_flattened_tree(next_position, keys, has_children, next_sibling);
      indentation -= 2;
      next_position = next_sibling[next_position];
    }
  }
}


int main(int argc, char ** argv) {
  int N = atoi(argv[1]);
  int Seed = atoi(argv[2]);
  srand(Seed);

  struct node_str * root = create_tree(N);

  print_tree(root);

  int * keys = (int *) malloc(sizeof(int) * N);
  int * has_children = (int *) malloc(sizeof(int) * N);
  int * next_sibling = (int *) malloc(sizeof(int) * N);

  flatten_tree(root, keys, has_children, next_sibling);
  print_flattened_tree(0, keys, has_children, next_sibling);

  cuda_compute(N, keys, has_children, next_sibling);

  return 0;
}